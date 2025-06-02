library(R6)
library(data.tree)
library(survival)
library(lpSolve)
library(osqp)

DipolarSurvivalTree_PolyKernel = R6Class(
  
  
  
  classname = "DipolarSurvivalTree_PolyKernel",
  
  public = list(
    
    # MAIN DATA
    traindata = NULL,
    originalrownames = character(),
    covariates = character(),
    time = NA_character_,
    censor = NA_character_,
    poms = matrix(NA_real_),
    pureweight = NA_real_,
    mixedweight = NA_real_,
    Kconstant = NA_real_,
    Kpoly_order = NA_integer_,
    K.bicovariates = matrix(NA_real_),
    
    # ALGORITHMS DATA
    quantiles = numeric(length = 2L), # cut-off quantiles to tag dipoles as "pure", "mixed", "neither"
    tolerance = NA_real_, # tolerance for reorientation algorithm convergence
    epsilon = NA_real_, # uniform margin for piecewise linear function
    kappa = NA_real_,
    ncovariatestosearch = NA_integer_,
    
    # TREE DATA
    nsize = NA_integer_, # terminal node size threshold
    
    # INITIALIZER
    initialize = function(
      traindata, time, censor, covariates,
      quantiles = c(.35,.65), tolerance = 10^-3, epsilon = 1, kappa = 1, nsize = 10L,
      pureweight = 1, mixedweight = 1,
      Kconstant = 1, Kpoly_order = 2,
      ncovariatestosearch = 10
    ) {
      # Don't allow empty data sets
      if (!is.data.frame(traindata) || nrow(traindata) == 0) {
        stop("\"traindata\" variable must be a nonempty data frame")
      }
      # Rename rows of "traindata" to a sequence of consecutive integers
      # Store the original row names in case they are required later
      self$originalrownames = rownames(traindata)
      rownames(traindata) = seq(1:nrow(traindata))
      # Make sure the survival and covariate variable names are actually valid
      for (name in c(time, censor, covariates)) {
        if (!(name %in% names(traindata))) {
          stop(
            paste(name, "is not a variable in the data", sep = " ")
          )
        }
      }
      # Bulk initiation of various variables
      self$traindata = traindata
      self$covariates = covariates
      self$time = time
      self$censor = censor
      self$quantiles = quantiles
      self$tolerance = tolerance
      self$epsilon = epsilon
      self$kappa = kappa
      self$nsize = nsize
      
      self$ncovariatestosearch = min(ncovariatestosearch, length(covariates))
      
      if (pureweight <= 0 || mixedweight <= 0) {
        stop("weights for dipoles must be positive")
      }
      self$pureweight = pureweight
      self$mixedweight = mixedweight
      # Create the "pure or mixed" matrix
      self$poms = self$createpoms()
      
      
      self$Kconstant = Kconstant
      self$Kpoly_order = Kpoly_order
      
      
      covariates.mat = as.matrix(traindata[, covariates])
      self$K.bicovariates = self$K.bi(covariates.mat)
    },
    
    
    K.bi = function(X) {
      (X%*%t(X) + self$Kconstant)^self$Kpoly_order
    },
    
    K = function(X, x) {
      (X%*%x + self$Kconstant)^self$Kpoly_order
    },
    
    # CREATE "PURE OR MIXED" MATRIX
    createpoms = function() {
      lowerQ = self$quantiles[1]
      upperQ = self$quantiles[2]
      # Make pairwise distance vector while accounting for censoring
      N = nrow(self$traindata)
      DT <- rep(NA_real_, ceiling(0.5*N*N))
      DTindex = 1 ## PREVIOUSLY HAD ERROR! R starts indices from 1
      for (i in 1:(N-1)) {
        # survival information of i-th data point
        xitime = self$traindata[i, self$time]
        xicensor = self$traindata[i, self$censor]
        for (j in (i+1):N) {
          # survival information of j-th data point
          xjtime = self$traindata[j, self$time]
          xjcensor = self$traindata[j, self$censor]
          # account for censoring as follows:
          if (
            (xicensor == 1) && (xjcensor == 1) ||
            ((xicensor == 0) && (xjcensor == 1) && (xitime > xjtime)) ||
            ((xicensor == 1) && (xjcensor == 0) && (xjtime > xitime))
          ) {
            DT[DTindex] = abs(xitime - xjtime)
            DTindex = DTindex + 1
          }
        }
      }
      # Get the lower and upper quantile cut-off thresholds
      QV = quantile(DT, c(lowerQ, upperQ), na.rm = TRUE)
      lQV = QV[[1]]
      uQV = QV[[2]]
      # Determines if a pair of covariates (xi, xj) is pure, mixed or neither
      # according to Kretowska 2017 and stores the results in a matrix
      poms = matrix(0,N,N)
      for (i in 1:(N-1)) {
        xitime = self$traindata[i, self$time]
        xicensor = self$traindata[i, self$censor]
        for (j in (i+1):N) {
          xjtime = self$traindata[j, self$time]
          xjcensor = self$traindata[j, self$censor]
          if (
            (xicensor == 1) && (xjcensor == 1) && (abs(xitime - xjtime) <= lQV)
          ) {
            poms[i, j] = 1
          } else if (
            ((xicensor == 1) && (xjcensor == 1) && (abs(xitime - xjtime) >= uQV)) ||
            ((xicensor == 0) && (xjcensor == 1) && ((xitime - xjtime) >= uQV)) ||
            ((xicensor == 1) && (xjcensor == 0) && ((xjtime - xitime) >= uQV))
          ) {
            poms[i, j] = -1
          }
        }
      }
      return(poms)
    },
    
    createtree = function(subset) {
      subX = self$traindata[subset, self$covariates]
      subtree = Node$new("Node0", data = subX)
      private$smallestweightpure = 1
      private$smallestweightmixed = 1
      self$createtreenode(subtree)
      subtree$Do(function(node) node$KMest <-self$KMcompute(node),
                 filterFun = isLeaf)
      subtree$Do(function(node) node$pureprop <-self$pureprop.compute(node))#,
      #           filterFun = isLeaf)
      return(subtree)
    },
    
    # CREATE OBLIQUE TREE RECURSIVELY
    createtreenode = function(node) {
      n = nrow(node$data)
      # Split only if node size is larger than threshold
      if (n > nsize) {
        # Instantiate book-keeping variables that are shared across many functions
        # These change at every iteration
        private$X = node$data
        # private$Z = as.matrix(cbind(1, private$X))
        private$subsetX = strtoi(rownames(private$X))
        private$K.XX = self$K.bicovariates[private$subsetX, private$subsetX]
        private$D = ncol(private$X)
        private$N = nrow(private$X)
        private$N2 = 2*private$N
        # Get the optimal quantities necessary for split
        w0_mupmdiff = self$optimal_w0_mupmdiff()
        node$opt_w0_mupmdiff <- w0_mupmdiff 
        lrnodes = self$lrnode.make(w0_mupmdiff)
        nl = nrow(lrnodes[[1]])
        nr = nrow(lrnodes[[2]])
        # Only grow the tree if both left and right nodes are nonempty
        if (nl > 0 && nr > 0) {
          splits <- lrnodes[[3]]
          lrstat <- self$lr.compute(splits)
          node$lrstat <- lrstat
          # label child nodes
          lcnodename = paste0(node$name,"l")
          rcnodename = paste0(node$name,"r")
          # grow the tree by adding child nodes to current node
          node$AddChild(lcnodename, data = lrnodes[[1]]) 
          node$AddChild(rcnodename, data = lrnodes[[2]])
          # recursively apply tree growing to each child
          self$createtreenode(node$children[[1]])
          self$createtreenode(node$children[[2]])
        }
      }
    },
    
    # GET OPTIMAL LAGRANGE VARIABLES OF CURRENT NODE DATA STORED IN private$X
    optimal_w0_mupmdiff = function() {
      private$subpoms = self$poms[private$subsetX, private$subsetX]
      w0_mupmdiff = private$optimizer()
      return(w0_mupmdiff)
    },
    
    # USE OPTIMAL LAGRANGE VARIABLES TO SPLIT CURRENT private$X
    lrnode.make = function(w0_mupmdiff){
      w0 = w0_mupmdiff$w0
      mupmdiff = w0_mupmdiff$mupmdiff
      splits = private$K.XX %*% mupmdiff + w0
      XL = private$X[splits < 0, ]
      XR = private$X[splits >= 0, ]
      return(list(XL, XR, splits))
    },
    
    # COMPUTE LOG RANK STATISTIC OF CURRENT private$X (WITH ITS SURVIVAL INFORMATION)
    lr.compute = function(splits){
      subdata = self$traindata[private$subsetX, ]
      Y = Surv(
        subdata[self$time][[1]], subdata[self$censor][[1]] == 1
      )
      lrstat = survdiff(Y ~ splits < 0)[[5]]
      return(lrstat)
    },
    
    
    # METHOD FOR ADDING KM ESTIMATES TO TERMINAL NODES 
    KMcompute = function(node) {
      X<-node$data 
      subsetX = strtoi(rownames(X))
      datasubset = self$traindata[subsetX,]
      Y<-Surv(datasubset[time][[1]],datasubset[censor][[1]]==1)
      kmfit=survfit(Y~1)
      return(kmfit) 
    },
    
    
    # PURE PROPORTION IN TERMINAL
    pureprop.compute = function(node) {
      subrows = strtoi(rownames(node$data))
      subpoms = self$poms[subrows, subrows]
      subtable = table(subpoms)
      pureprop = if(is.na(subtable["1"]) && is.na(subtable["-1"])) {
        NA
      } else if(is.na(subtable["-1"])) {
        1
      } else {
        subtable["1"]/(subtable["1"] + subtable["-1"])
      }
      return(pureprop)
    },
    
    
    # ADD PRUNING CRITERIA NUMBERS
    pruningstat = function(node) {
      GTh<-sum(node$Get("lrstat",filterFun=isNotLeaf))
      Sh <- node$totalCount - node$leafCount
      gh <- GTh/Sh
    },
    
    
    
    # PREDICT TIME
    predicttime = function(testdata, trainedtree) {
      testdataX = testdata[self$covariates]
      predictedtimes = rep(NA, nrow(testdataX))
      testdataX.mat = as.matrix(testdataX)
      for (i in 1:nrow(testdataX.mat)) {
        private$testx = testdataX.mat[i,]
        self$predrec(trainedtree)
        predictedtimes[i] = private$predtesttime
      }
      return(predictedtimes)
    },
    
    
    # RECURSIVELY TRAVERSE FOR PREDICTED TIME
    predrec = function(node) {
      if (isLeaf(node)) {
        km = node$KMest
        kmmed = summary(km)$table["median"][[1]]
        pred = if (is.na(kmmed)) {
          #summary(km)$table["*rmean"][[1]]
          mean(self$traindata[rownames(node$data),self$time])
        } else {
          kmmed
        }
        private$predtesttime = pred
      } else {
        X = as.matrix(node$data)
        w0 = node$opt_w0_mupmdiff$w0
        mupmdiff = node$opt_w0_mupmdiff$mupmdiff
        split = (t(mupmdiff) %*% self$K(X, private$testx)) + w0
        if (split < 0) {
          self$predrec(node$children[[1]])
        } else {
          self$predrec(node$children[[2]])
        }
      }
    },
    
    
    
    
    # CONCORDANCE INDEX
    cindex = function(actualtime, predictedtime, censor) {
      
      N = length(actualtime)
      matches = 0
      legitcomparisons = 0
      for(i in 1:N) {
        for(j in 1:N) {
          actualless = (actualtime[i] < actualtime[j])
          predictedless = (predictedtime[i] < predictedtime[j])
          concorded = actualless && predictedless
          
          notsamenode = predictedtime[i] != predictedtime[j]
          
          matches = matches+censor[i]*as.numeric(
            concorded && notsamenode
          )
          
          legitcomparisons = legitcomparisons+censor[i]*as.numeric(
            actualless && notsamenode
          )
        }
      }
      matches/legitcomparisons
    },
    
    
    ridgedipolarcriterion = function(betapm, K.XX, mupmdiff, w0) {
      wX_plus_w0_kernelized = t(mupmdiff) %*% K.XX  + w0
      
      ww_kernelized = t(mupmdiff) %*% K.XX %*% mupmdiff
      
      phipm_kernelized = c(
        pmax(0, self$epsilon - wX_plus_w0_kernelized),
        pmax(0, self$epsilon + wX_plus_w0_kernelized)
      )
      
      Psi_betapm_w_w0_kernelized = sum(betapm * phipm_kernelized)
      
      return(0.5 * ww_kernelized + self$kappa * Psi_betapm_w_w0_kernelized)
    },
    
    
    
    dipolarcriterion = function(betapm, X, w, w0) {
      wX_plus_w0 = (X %*% w)  + w0
      
      ww = sum(w*w)
      
      phipm = c(
        pmax(0, self$epsilon - wX_plus_w0),
        pmax(0, self$epsilon + wX_plus_w0)
      )
      
      Psi_betapm_w_w0 = sum(betapm * phipm)
      
      return(0.5 * ww + self$kappa * Psi_betapm_w_w0)
    }
    
    
  ),
  
  
  
  
  
  private = list(
    
    # # CONSTANT VARIABLE FOR INTERNAL USE ONLY
    # vinitial.default.index = NA_integer_,
    
    # SHARED BOOK-KEEPING VARIABLES FOR INTERNAL USE ONLY
    # THESE CHANGE FREQUENTLY AND UNPREDICTABLY
    # Z = NULL,
    X = NULL,
    subsetX = integer(),
    K.XX = matrix(NA_real_),
    subpoms = matrix(NA_real_),
    
    D = NA_integer_,
    N = NA_integer_,
    N2 = NA_integer_,
    smallestweightpure = NA_real_,
    smallestweightmixed = NA_real_,
    
    testx = numeric(),
    predtesttime = NA_real_,
    
    # COEFFICIENT CALCULATOR FOR THE DIPOLAR CRITERION FUNCTION
    betapm = function(Zv_kernelized){
      phip <- rep(0, private$N)
      phim <- rep(0, private$N)
      for (i in 1:(private$N - 1)) {
        for (j in (i+1):private$N) {
          if (private$subpoms[i, j] == 1) {
            # if we EXPECT both zi and zj on +ve side of v:
            if (Zv_kernelized[i] + Zv_kernelized[j] > 0) {
              phip[i] = phip[i] + private$smallestweightpure
              phip[j] = phip[j] + private$smallestweightpure
              # otherwise:
            } else {
              phim[i] = phim[i] + private$smallestweightpure
              phim[j] = phim[j] + private$smallestweightpure
            }
          } else if (private$subpoms[i, j] == -1) {
            # if we EXPECT zi on +ve side of v
            # while zj on -ve side of v:
            if (Zv_kernelized[i] - Zv_kernelized[j] > 0) {
              phip[i] = phip[i] + private$smallestweightmixed
              phim[j] = phim[j] + private$smallestweightmixed
              # otherwise:
            } else {
              phim[i] = phim[i] + private$smallestweightmixed
              phip[j] = phip[j] + private$smallestweightmixed
            }
          }
        }
      }
      return(
        list(Zv_kernelized, c(phip, phim))
      )
    },
    
    # w0.manual = function(mupm, betapm, mupmdiff) {
    #   N = length(mupmdiff)
    #   
    #   gammapm = mupm/betapm
    #   
    #   betapm.non0index = which(betapm != 0)
    #   
    #   closenesserror = mean(gammapm[betapm.non0index])/1000 # tolerance to determine "close enough" to 0 or kappa
    #   
    #   gammapm.inactive = sapply(betapm.non0index, function(i) {
    #     gammapm.i = gammapm[i]
    #     # because of numerical imprecision, some gammapm's that should be on the 
    #     # constraint boundaries could be slightly off of them; this is a brute
    #     # force way to filter those out using an absolute tolerance
    #     # gammapm.i > 0 && gammapm.i < kappa && abs(gammapm.i) > closenesserror && abs(kappa - gammapm.i) > tol
    #     gammapm.i > closenesserror && (self$kappa - gammapm.i) > closenesserror
    #   })
    #   gammapm.inactiveindex = betapm.non0index[gammapm.inactive]
    #   
    #   w0s = sapply(gammapm.inactiveindex, function(k) {
    #     if(k <= N) {
    #       return(self$epsilon - sum(mupmdiff*private$K.XX[k,]))
    #     } else {
    #       k = k %% N
    #       if(k == 0) {k = N}
    #       return(-self$epsilon - sum(mupmdiff*private$K.XX[k,]))
    #     }
    #   })
    #   
    #   if(sum(gammapm.inactive) == 0) {
    #     print("WARNING! w0 ISSUE!")
    #     return(0)
    #   } else{
    #     return(mean(w0s))
    #   }
    # },
    
    
    # OPTIMIZER TO FIND OPTIMAL SPLITTING VECTOR
    optimizer = function() {
      settings <- osqpSettings(verbose = FALSE, max_iter = 10000L,
                               eps_abs = 1e-6,
                               eps_rel = 1e-6,
                               adaptive_rho_interval = 25,
                               warm_start = TRUE)
      
      epsvec = rep(self$epsilon, private$N)
      q = -c(epsvec, epsvec) # for max negate coefficients
      P = kronecker(
        matrix(c(1, -1, -1, 1), byrow = FALSE, nrow = 2),
        private$K.XX
      ) # for max negate coefficients
      P = Matrix::Matrix(P, sparse = TRUE)
      
      Aeq = t(rep(c(1, -1), c(private$N, private$N)))
      leq = ueq = 0
      
      A = rbind(Aeq, diag(private$N2))
      A = Matrix::Matrix(A, sparse = TRUE)
      l = c(leq, rep(0, private$N2))
      
      model = osqp(P = P, q = q, A = A, l = l, u = NULL, pars=settings)
      
      covariatestosearch = sample(1:private$D, 
                                  size = self$ncovariatestosearch,
                                  replace = FALSE)
      ncovariatestosearch3 = 3*self$ncovariatestosearch
      dualsols = numeric(ncovariatestosearch3)
      ress = vector("list", length = ncovariatestosearch3)
      # try out ncovariatestosearch random covariates
      for (j.index in 1:self$ncovariatestosearch) {
        j = covariatestosearch[j.index]
        w0s = -quantile(private$X[, j], probs = c(.25, .5, .75))
        w = rep(0, private$D)
        w[j] = 1
        
        for(jj.index in 1:3) {
          Zv <- as.matrix(private$X) %*% w + w0s[jj.index]
          betapm = private$betapm(Zv)[[2]]
          
          u = c(ueq, self$kappa*betapm)
          model$Update(u = u)
          res = model$Solve()
          
          jjj.index = 3*(j.index - 1) + jj.index
          dualsols[jjj.index] = -res$info$obj_val #; print(j); print(dualsols[j])
          ress[[jjj.index]] = res
        }
      }
      # whichever produced the smallest objective after optimizing ONCE,
      # choose that to continue recursive orientation
      jjj.index.min = which.min(dualsols)
      
      # vinitial = private$vinitial.default()
      # phipm.list <- private$phipmcount(vinitial)
      # betapm = phipm.list[[2]]
      
      res = ress[[jjj.index.min]]
      
      mupm = res$x
      mup = mupm[1:private$N]; mum = mupm[(private$N + 1):(private$N2)]
      mupmdiff = mup - mum
      
      #w = c(t(mup - mum) %*% X)
      
      dualsol = -res$info$obj_val # minus because we had to invert signs for max earlier
      
      #w0 = private$w0(mupm, betapm, mupmdiff)
      
      w0 = res$y[1]
      
      #print(w0)
      #print(res$y[1])
      
      i = 0
      
      dualsol.old = dualsol
      
      dualsol.min = dualsol
      mupmdiff.min = mupmdiff
      w0.min = w0
      
      # print(
      #   paste("dual:", dualsol)
      # )
      # print(
      #   paste("primal a opt:", self$ridgedipolarcriterion(betapm, private$K.XX, mupmdiff, w0))
      # )
      repeat {
        Zv_kernelized <- private$K.XX %*% mupmdiff + w0
        betapm = private$betapm(Zv_kernelized)[[2]]
        #print(betapm)
        # print(
        #   paste("primal a opt a reorient:", self$ridgedipolarcriterion(betapm, private$K.XX, mupmdiff, w0))
        # )
        
        u = c(ueq, self$kappa*betapm)
        model$Update(u = u)
        res = model$Solve()
        
        mupm = res$x
        mup = mupm[1:private$N]; mum = mupm[(private$N + 1):(private$N2)]
        mupmdiff = mup - mum
        
        dualsol = -res$info$obj_val # minus because we had to invert signs for max earlier
        
        w0 = res$y[1]
        
        if(dualsol <= dualsol.min) {
          dualsol.min = dualsol
          mupmdiff.min = mupmdiff
          w0.min = w0
        }
        
        
        # print(
        #   paste("dual:", dualsol)
        # )
        # print(
        #   paste("primal a opt:", self$ridgedipolarcriterion(betapm, private$K.XX, mupmdiff, w0))
        # )
        i = i + 1
        if(abs(dualsol - dualsol.old) < max(self$tolerance * dualsol, 0.1*self$tolerance)) {
          # print("END OF REORIENT")
          # print(" ")
          #print(dualsol.min)
          #print(dualsol)
          break
        }
        
        if(i > 250) {
          print("max iter 250 reached")
          #print(dualsol.min)
          break
        }
        dualsol.old = dualsol
      }
      
      return(
        list(w0 = w0.min, mupmdiff = mupmdiff.min)
      )
    }
    
  )
  
)


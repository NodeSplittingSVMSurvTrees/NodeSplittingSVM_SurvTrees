intvarfun = function(data,covariates)
{
  int.mean = colMeans(alldata[covariates])  
  int.var =mean(rowSums((alldata[covariates] -int.mean)^2))
  return(int.var)
}
choosekappa <- function(alldata.train,etalist=c(),model){
  CIMs = c(); BScores=c();j=0; 
  #Find index for k folds for data set 
  k=5; 
  n=nrow(alldata.train)
  a = rep(floor(n/k),k) 
  b = rep(0:1,times=c(k-n%%k,n%%k)) 
  c = a + b 
  folds<-split(1:n, sample(rep(1:k,times=c)))
  for (eta in etalist){
    j=j+1; 
    if (model == "l") {
      Dipolar.model <- DipolarSurvivalTree_PolyKernel$new(
        alldata.train, time, censor, covariates,
        quantiles, tolerance, epsilon, kappa=exp(eta), nsize,
        pureweight=1, mixedweight=1,Kconstant=0,Kpoly_order =1
      )
    } else if (model == "q") {
      Dipolar.model <- DipolarSurvivalTree_PolyKernel$new(
        alldata, time, censor, covariates,
        quantiles, tolerance, epsilon, kappa=exp(eta), nsize,
        pureweight=1, mixedweight=1,Kconstant=1,Kpoly_order =2
      )
    } else {
      ksigma=intvarfun(alldata.train,covariates)
      Dipolar.model <- DipolarSurvivalTree_GaussKernel$new(
        alldata.train, time, censor, covariates,
        quantiles, tolerance, epsilon, kappa=exp(eta), nsize,
        pureweight=1, mixedweight=1, Ksigma = ksigma
      )
    }
  CIndex = c(); BScore = c(); 
  for (i in 1:5){
    trainsubset<-setdiff(1:nrow(alldata.train),folds[[i]])
    dipolartree.kfold<-Dipolar.model$createtree(trainsubset)
    dipolartree.kfold.prune<-bootstrapPruning(dipolartree.kfold,Dipolar.model,2.8)[[5]]
    alldatatest=Dipolar.model$traindata[folds[[i]],]
    predictedtime2 = Dipolar.model$predicttime(alldatatest,dipolartree.kfold.prune)
    actualtime = alldatatest$stop
    #actualtime = alldatatest$survt
    CIndex[i] = Dipolar.model$cindex(actualtime,predictedtime2,alldatatest$status)
    IBSrange = gmsfun(alldatatest$stop,alldatatest$status,2)
    BScore[i]=Int_Brierscore(folds[[i]],dipolartree.kfold.prune,Dipolar.model,IBSrange,0)
  }
  CIMs[j] = mean(CIndex)
  BScores[j] = mean(BScore)
  }

results=list(CIMs,BScores,etalist[which.max(CIMs)])
return(results)
} 

#FUNCTIONS FOR SIMULATIONS ----
#Split Complexity 
splitcomplexity<- function(node,alpha=2.8) {
  if (node$totalCount==1) {GTh = node$lrstat} 
  else {GTh<-sum(node$Get("lrstat",filterFun=isNotLeaf))}
  Sh <- node$totalCount - node$leafCount
  gh <- GTh-alpha*Sh
}

#Function for IBSRange 
gmsfun = function(cforgms,status,k){
  cforgms[status==1]
  cforgms=sort(cforgms[!duplicated(cforgms)])
  n=length(cforgms); a=c()
  for (i in 1:(n-1)){
    a[i]=sqrt(cforgms[i]*cforgms[i+1])
  }
  a=a[seq(1,length(a),k)]
  return(a)
}

#Make k fold partition
foldsfun =function(data,k){
  n=nrow(data)
  a = rep(floor(n/k),k) 
  b = rep(0:1,times=c(k-n%%k,n%%k)) 
  c = a + b 
  folds<-split(1:n, sample(rep(1:k,times=c)))
  return(folds)
}



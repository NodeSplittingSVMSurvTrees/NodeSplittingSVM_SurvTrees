predfunction = function(Dipolar.model, ourtree,x,ptimes,splitrule){
  #splitrule==0 means use Z %*% v < 0 for splitting rule, splitrule==1 means use 
  # Z %*% v <= 0 for splitting rule. 
  if (splitrule == 0) {
    node = ourtree
    while (isNotLeaf(node)){
      X = as.matrix(node$data)
      w0 = node$opt_w0_mupmdiff$w0
      mupmdiff = node$opt_w0_mupmdiff$mupmdiff
      splittest = ((t(mupmdiff) %*% Dipolar.model$K(X, x)) + w0 < 0)
      if (splittest){
        node=node$children[[1]]
      } else {
        node=node$children[[2]]
      } 
    }
  } else {
    z = c(1,x)
    node = ourtree
    while (isNotLeaf(node)){
      v = node$optv
      splittest = (v%*%z<=0)
      if (splittest){
        node=node$children[[1]]
      } else {
        node=node$children[[2]]
      } 
    }
  }
  
  outnode=node$name
  pred=summary(node$KMest,times=ptimes,extend=TRUE)[[6]]
  return(list(outnode,pred))
}


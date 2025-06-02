Brierscore2 = function(fold,ourtree,ourmodel,t,splitrule){
testdata=ourmodel$traindata[fold,]
n = nrow(testdata)
kmCfit = survfit(Surv(testdata$stop,testdata$status==0)~1)
#kmCfit = survfit(Surv(ourmodel$traindata$stop,ourmodel$traindata$status==0)~1)
X = testdata[ourmodel$covariates]
s=rep(0,n)
for (i in 1:n){
  pred = predfunction(ourmodel,ourtree,as.numeric(as.vector(X[i,])),t,splitrule)[[2]]
  T = testdata['stop'][i,]; d = testdata['status'][i,]
  GT = summary(kmCfit,times=T)[[6]]; 
  if (GT==0){
    y=summary(kmCfit)[[6]]
    if (!all(y==0)){
      GT = min(y[y!=min(y)])
      }
    else GT=1
    }
  Gt = summary(kmCfit,times=t,extend=TRUE)[[6]]
  s[i] = pred^2*1*(T<=t)*d/GT +  (1-pred)^2*1*(T>t)/Gt
  
}
s=sum(s)/n; 
}


trapz = function(x, y) 
{ # computes the integral of y with respect to x using trapezoidal integration. 
  idx = 2:length(x)
  return (as.double( (x[idx] - x[idx-1]) %*% (y[idx] + y[idx-1])) / 2)
}

Int_Brierscore = function(fold,ourtree,ourmodel,ts,splitrule){
  a = length(ts); BSt = rep(0,a)
  for (i in 1:a){
    BSt[i]=Brierscore2(fold,ourtree,ourmodel,ts[i],splitrule)
  } 
IBS = trapz(ts,BSt)/max(ts)
}


bootstrapPruning <- function(intree,dipolarmodel,alpha) {

#bts[1] = number of bootstrap samples, bts[2] = percentage size of bootstrap relative to sample
#splitrule==0 means use Z %*% v < 0 for splitting rule, splitrule==1 means use 
# Z %*% v <= 0 for splitting rule. 

#I) Get list of pruned trees from all data (T1 > T2 > ... > Tm)  tree and list of pruning stats (alpha_1 < ... < alpha_m)

#Add pruning criteria numbers function 
pruningstat<- function(node) {
  GTh<-sum(node$Get("lrstat",filterFun=isNotLeaf))
  Sh <- node$totalCount - node$leafCount
  gh <- GTh/Sh
}
alldata=dipolarmodel$traindata[strtoi(rownames(intree$data)),]
#function for pruning all children from node 
prunenodesfun = function(opttree, prunenode){
  prunenodel = paste0(prunenode,'l')
  Prune(opttree, pruneFun = function(x) x$name!=prunenodel)
  prunenoder = paste0(prunenode,'r')
  Prune(opttree, pruneFun = function(x) x$name!=prunenoder)
}


#function for getting list of pruned trees and pruning stats while pruning off 
#sequence of weakest links 

lrpruningfun = function(treetoprune) {
  pstats = c()
  opttree = list()
  opttree[[1]] = Clone(treetoprune)
  opttree[[1]]$Do(function(node) node$prnstat<- pruningstat(node), filterFun = isNotLeaf)  
  pstats[1] = min(opttree[[1]]$Get("prnstat", filterFun = isNotLeaf)) 
  i=0
  while (opttree[[i+1]]$totalCount != 1){
    i = i+1 
    prunenode=names(which.min(opttree[[i]]$Get("prnstat", filterFun = isNotLeaf)))
    prunenodesfun(opttree[[i]],prunenode)
    opttree[[i+1]] = Clone(opttree[[i]])
    opttree[[i+1]]$Do(function(node) node$prnstat<- pruningstat(node), filterFun = isNotLeaf) 
    if (opttree[[i]]$totalCount != 1){
      pstats[i+1] = min(opttree[[i+1]]$Get("prnstat", filterFun = isNotLeaf))
    }
  }
  opttree=opttree[-(i+1)]
  a = list(pstats,opttree)
  return(a)
} 
#Account for the case where the input tree has a single node 
if (intree$totalCount==1) {
  results = list("No pruning, input tree is single node","No pruning, input tree is single node",1,intree)
  return(results)
}
else {
listprune = lrpruningfun(intree)
listprunestat = listprune[[1]]
listprunetree = c(intree,listprune[[2]])
a = listprunestat>alpha
if (sum(a)==0) {index=length(a)+1} else 
  {index=min(which(a))}

KMcompute = function(node) {
  X<-node$data
  subsetX = rownames(X)
  datasubset = alldata[subsetX,]
  Y<-Surv(datasubset[time][[1]],datasubset[censor][[1]]==1)
  kmfit=survfit(Y~1)
  return(kmfit)
}
finaltree=listprunetree[index][[1]]
finaltree$Do(function(node) node$KMest=KMcompute(node), filterFun = isLeaf)
results = list(listprunestat,listprunetree,index,finaltree)
return(results)
}
}



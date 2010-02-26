util.sql.groupby <-
function(data, groupCols, sep="-", varCol, type=c("max", "min", "sum", "mean"), na.rm=FALSE){

  m = match(c("max", "min", "sum", "mean"), type, 0)
  matchFun = which(m==1)
  
  
  if(length(groupCols)>1){
    key = util.matrix.cat(data, cols=groupCols, sep=sep)
  }else{
    key = data[,groupCols]
  }

  keyVal = unique(key)
  grpMax=NULL
  grpMin=NULL
  grpSum=NULL
  grpMean= NULL
  varVec=NULL
  for( i in keyVal){
    filter = key==i
    varVec = unlist(data[filter, varCol])
    for(j in matchFun){
      if(j==1) grpMax = c(grpMax, max(varVec, na.rm=na.rm))
      if(j==2) grpMin = c(grpMin, min(varVec, na.rm=na.rm))
      if(j==3) grpSum = c(grpSum, sum(varVec, na.rm=na.rm))
      if(j==4) grpMean = c(grpMean, mean(varVec, na.rm=na.rm))
    }

  }
  
  re = NULL
  re = cbind(re, key=keyVal) 
  for (j in matchFun){
    if(j==1) re = cbind(re, max=grpMax)
    if(j==2) re = cbind(re, min=grpMin)
    if(j==3) re = cbind(re, sum=grpSum)
    if(j==4) re = cbind(re, mean=grpMean)
  }

  return(re)
}


util.matrix.rmSparseRow <-
function(ma, colChecked, minNumInCol=1){

  ifD = F
  checked = ma[,colChecked]
  if(ifD) print(checked)

  len = length(checked)
  rmList = NULL
  reList = NULL
  for( i in 1:len ){
    if(ifD) print(checked[i])
    if( checked[i] <= minNumInCol ){
      rmList = c(rmList, i)
    }else{
      reList = c(reList, i)

    }
  }
  if(ifD) print(rmList)
  if(is.null(rmList)) return(NULL)
  len = length(rmList)
  re = ma
    for( j in 1:len ){
      re = t(util.matrix.delCol (t(re), rmList[j]))
      if(ifD) print(re)
      rmList = rmList - 1
      if(ifD) print(rmList)
    }

  return(list(re, reList))
}


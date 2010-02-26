util.matrix.col.shuffle <-
function(ma){

  if(!is.null(dim(ma)))  return(ma)

  colNum = dim(ma)[2]
 
  filterSeq = matrix(1:colNum, ncol=2)
  filterSeq = as.vector(t(filterSeq))

  re = ma[,filterSeq]

  return(re)
}


util.matrix.col.shuffle2 <-
function(ma1, ma2){

  re = cbind(ma1, ma2)

  colNum = 1
  if(!is.null(dim(ma1))){
    colNum = dim(ma1)[2]
  }
  filterSeq = matrix(1:(2*colNum), ncol=2)
  filterSeq = as.vector(t(filterSeq))

  re = re[,filterSeq]

  return(re)
}


util.matrix.clone <-
function(ma, n, rowAppend=T){
  rnm = dim(ma)[1]
  cnm = dim(ma)[2]
  if(rowAppend){
    t1 = replicate(n, t(ma))
    re = t(matrix(t1, nrow=cnm, byrow=F))
    re
  }else{
    t1 = replicate(n, ma)
    re = matrix(t1, nrow=rnm, byrow=F)
    re
  }
  return(re)
}


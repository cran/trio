exhaustHapExp <-
function(lociCt=3, re.str=T,  snpCoding=c(1,2)){

  ma = matrix(snpCoding, ncol=1)

  ## grow the matrix on both sides
  if(lociCt==1) return(list(hap=ma, hapStr=as.character(ma)))
  for( i in 1:(lociCt-1) ){
    growing = util.matrix.clone(ma, 2)    
    addedBit = c(rep(snpCoding[1], times=2^i), rep(snpCoding[2], times=2^i))
    ma = cbind(growing, " "=addedBit)
  }
  if(re.str) {
    hapStr = util.matrix.cat(ma, 1:lociCt, sep="")
    return(list(hap=ma, hapStr=hapStr))
  }else{
    return(list(hap=ma))
  }
  
}


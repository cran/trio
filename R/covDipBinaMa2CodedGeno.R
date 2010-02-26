covDipBinaMa2CodedGeno <-
function(dip1, dip2, subjectCt, snpCoding=c(0,1,2,3), snpBase=c(0,1,2)){

  dipSum = dip1+dip2

  ##  print(dipSum)
  a1 = sum(snpBase[2:3] * c(2,0)) # 11-> to 1
  a2 = sum(snpBase[2:3] * c(1,1)) # 12-> to 3
  a3 = sum(snpBase[2:3] * c(0,2)) # 22-> to 2

  snp1d.f = factor(dipSum, levels=c(a1, a3, a2), labels = snpCoding[2:4])
  snp1d.f = as.character(snp1d.f)

  if(is.null(dim(dipSum))){
    ncol = length(dipSum)
    nrow = 1
  }else{
    ncol = dim(dipSum)[2]
    nrow = dim(dipSum)[1]
  }

  codedGeno = matrix(as.integer(snp1d.f), nrow=subjectCt, ncol=ncol)

  
##   codedGeno = matrix(NA, nrow=subjectCt, ncol=ncol)
##   for( i in 1:dim(dipSum)[1]){
##     for( j in 1:dim(dipSum)[2]){
##        if (dipSum[i,j] == 2*snpBase[1]) codedGeno[i,j] = snpCoding[1]
##        if (dipSum[i,j] == 2*snpBase[2]) codedGeno[i,j] = snpCoding[2]
##        if (dipSum[i,j] == 2*snpBase[3]) codedGeno[i,j] = snpCoding[3]
##        if (dipSum[i,j] == snpBase[2]+snpBase[3]) codedGeno[i,j] = snpCoding[4]
##     }
##   }
  return(codedGeno)
}


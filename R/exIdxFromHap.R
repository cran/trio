exIdxFromHap <-
function(lociCt){
  ## know the internal sequency of haplotype
  digitMap1 = matrix(NA, ncol=lociCt, nrow=2^(lociCt-1) );
  digitMap2 = matrix(NA, ncol=lociCt, nrow=2^(lociCt-1) );

  for( n in 1:lociCt ){
     col = seq.int(from=1, to=2^lociCt, by = 2^n)
     col2 = col+2^(n-1)

     colNext = col
     colNext2 = col2

     for( j in 1:(2^(n-1))){
       if( j > 1){
         newCol = col + j - 1
         colNext = rbind(colNext, newCol)

         newCol2 = col2 + j - 1
         colNext2 = rbind(colNext2, newCol2)
         ## print(colNext)
       }
     }

     digitMap1[,n]=as.vector(colNext)
     digitMap2[,n]=as.vector(colNext2)
   }

  
  return (list(digitMap1=digitMap1, digitMap2 = digitMap2))
}


findMissing <-
function(df, is.1digit=F, snpStartLeftIndex, snpEndRightIndex, dig1Code=c(0, 1, 2, 3), dig2Code=c(0, 1, 2) ){
  ## use c(0,1,3,2) as the 1-digit coding, if not, exchange it
  snpEndLeftIndex = ncol(df)-snpEndRightIndex+1

  if(!is.1digit){
     snpNum = (snpEndLeftIndex - snpStartLeftIndex +1) /2
     snps = df[, snpStartLeftIndex:snpEndLeftIndex]

     snp1digit = exchangeDigit(ma=snps, cols=c(1,snpNum*2), dig1Code=c(0, 1, 3, 2), dig2Code =dig2Code, action=c("2to1"))
 
   }else{
     snpNum = snpEndLeftIndex - snpStartLeftIndex +1
     snp1digit = df[, snpStartLeftIndex:snpEndLeftIndex]
     if(min(dig1Code==c(0,1,3,2))==0){
       snp1digit = apply(snp1digit, 1:2,  FUN= util.vec.replace, orignal = dig1Code, replaceBy=c(0,1,3,2))
     }
   }
   nrow = nrow(df)
   missingPos = which(snp1digit==0)

   if(length(missingPos)==0) return(NULL)

   missingCol = (missingPos - (missingPos %% nrow))/nrow + 1

   miss.cord = matrix(NA, nrow=length(missingPos), ncol=2)

   miss.cord[,1]= (missingPos %% nrow)
   ## adjust the one at the last row
   miss.cord[  miss.cord[,1]==0, 1]= nrow
  
   miss.cord[,2]= (missingPos - miss.cord[,1])/nrow + 1

   colnames(miss.cord)=c("row", "snp")
   return(miss.cord)

}


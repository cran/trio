util.it.triMatch2 <-
function(dipIdx, len, re.homo=F){
       ## rely on the structure of tri idx pair
       ## rely on the order to the pair

       upper = .5*len*(len-1)
       if ((dipIdx)>upper){
         if(re.homo){
           return(c(rep(dipIdx-upper, 2), 0))
         }else{
           return(rep(dipIdx-upper, 2))
         }
       }

       cutoff = (len-1):1
       cutoff = cumsum(cutoff)

       block = qing.cut(dipIdx, cutPt=cutoff, cutPt.ordered = T, right.include=T)
       if(re.homo){
         return(c(block, dipIdx - c(0, cutoff)[block] + block, 1))
       }else{
         return(c(block, dipIdx - c(0, cutoff)[block] + block))
       }
#       if (block==1){
#         return( c(block, dipIdx+1))
#       }else{
#         return( c(block, dipIdx - cutoff[block-1] + block))
#       }

}


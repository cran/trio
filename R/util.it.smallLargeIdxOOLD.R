util.it.smallLargeIdxOOLD <-
function(idx, keep.same=F){
    ct = length(idx)
    idxPair=NULL
    if(keep.same){
     for ( i in 1:ct ){
       for ( j in 1:ct){
         #if(i<j) IdxPair = rbind(hapIdxCorner, c(hap1=i, hap2=j))
         if(i<=j) idxPair = rbind(idxPair, c(idx1=i, idx2=j))
       }
     }
    return(idxPair)
   }else{
     for ( i in 1:ct ){
       for ( j in 1:ct){
         #if(i<j) IdxPair = rbind(hapIdxCorner, c(hap1=i, hap2=j))
         if(i<j) idxPair = rbind(idxPair, c(idx1=i, idx2=j))
       }
     }
     return(idxPair)
   }
}


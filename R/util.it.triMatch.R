util.it.triMatch <-
function(idxPair, len){
  ## rely on the structure of tri idx pair
  ## rely on the order to the pair
  if(idxPair[1]-idxPair[2]==0){
     ## two idx are the same
     endIdx = (len^2-len)/2
              rowIdx = endIdx + idxPair[1]
   }else{
              # two idx are not the same, old approach
              #idxLen = c(0, (len-1):1)
              #idxLen.upper = idxLen[1:idxPair[1]]
              #idx.added = idxPair[2]-idxPair[1]
              #rowIdx = sum(idxLen.upper)+idx.added

              # new math approach
              idx.added = idxPair[2]-idxPair[1]
              rowIdx = len/2*(len-1)-(1+len-idxPair[1])/2*(len-idxPair[1]) + idx.added
   }
   return(rowIdx)
}


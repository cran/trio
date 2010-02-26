filterHapIdx2SuperHapSet <-
function( bkIdxes, hapIdxes, hapCts){

  it = lapply(hapCts, 1, FUN=seq, from=1)
  # first build the whole list
  idxComb = qExpandTable(listOfFactor = it, removedRowIdx=NULL, re.row=F )
  
  # eliminate the impossible rows.
  left = 1: (cumprod(hapCts)[length(hapCts)]  )
  for ( i in 1:length(bkIdxes)){
    matchedVal = hapIdxes[[i]]
    set = left[ is.element(idxComb [left ,bkIdxes[i]], matchedVal) ]
    #print(set)
    left = set
  }
  return(set)

}


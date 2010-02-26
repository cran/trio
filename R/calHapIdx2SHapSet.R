calHapIdx2SHapSet <-
function( bkIdx, hapIdx, hapCts){

  ifD = F
  
  ## find the number of row for each stratum
  cumRows = c(1, hapCts[1:bkIdx])
  cumRows = cumprod(cumRows)

  it = cumRows[bkIdx]
  matchSet = 1:it
  if(ifD) print(paste("matchSet=", paste(matchSet, collapse=";")))
  set = it*(hapIdx-1)+matchSet
  if(ifD) print(paste("set=", paste(set, collapse=";")))
  
  strataOffset = seq.int(from=0, to=cumprod(hapCts)[length(hapCts)]-1, by = it*hapCts[bkIdx])
  if(ifD) print(paste("strataOffset=", paste(strataOffset, collapse=";")))
  set = rep(strataOffset, each=it) + set;
  return(set)

}


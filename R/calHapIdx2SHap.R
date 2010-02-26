calHapIdx2SHap <-
function(hapIdxes, hapCts){

  # return one index
  bkCt = length(hapCts)

  hapCts= c(1, hapCts[-bkCt])
  cumRows = cumprod(hapCts)
  enuStart = cumRows*(hapIdxes-1)
  idx = cumsum(enuStart)[bkCt] + 1

  return(idx)

}


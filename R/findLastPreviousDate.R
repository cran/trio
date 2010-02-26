findLastPreviousDate <-
function(timeVec, probVec, cutoff){
  if(max(timeVec)< cutoff) return(NA)
  filter = which(timeVec<= cutoff)
  re = probVec[filter[length(filter)]]
  return(re)
}


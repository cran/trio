util.matrix.colComp <-
function(ma, values, operator="and"){
  ifD = F
  vCt = length(values)

  ## if "and" operator is chosen, the original list will decrease to speed up the comparison
  matchId = 1:(dim(ma)[1])
  if(ifD) print(matchId)
  for (i in 1:vCt ){
    if(ifD) print(i)
    compIdx = which(ma[matchId,i]==values[i])
    matchId = matchId[compIdx]
    if(ifD) print(paste("compIdx=", compIdx, collapse=", ", sep=""))
    if(ifD) print(paste("matchId=", matchId, collapse=", ", sep=""))
    
    if(length(matchId)==0) return (NULL)
  }
  return(matchId)

  ## if "or" operator is chose, the original list will remain the same to capture all that meet it.
  ## not implemented
}


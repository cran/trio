util.vec.matchVec <-
function(vec, benchMark, vecLen, benchLen){

  dup = vecLen/benchLen

  bench = rep(benchMark, dup)

  matchset = vec == bench

  re = NULL
  for( i in 0:(dup-1)){
    pos = i*benchLen + 1
    re = c(re, as.logical(min(matchset[pos:(pos+benchLen-1)])))
  }

  return(re)
}


util.vec.matchVecIdx <-
function(vec, benchMark, vecLen, benchLen){

  aa = util.vec.matchVec(vec, benchMark, vecLen, benchLen)
 
  return(which(aa))
}


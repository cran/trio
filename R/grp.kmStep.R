grp.kmStep <-
function(timeVec, probVec, lty=1, lwd=1,  col ="black", plot.dot=F){

  len = length(timeVec)

  for( i in 1:(len-1)){
    segments(timeVec[i], probVec[i], timeVec[i+1], probVec[i], col=col, lty=lty, lwd=lwd)
    segments(timeVec[i+1], probVec[i], timeVec[i+1], probVec[i+1], col=col, lty=lty, lwd=lwd)
    if(plot.dot) points(timeVec[i+1], probVec[i], pch=1, col="black")
  }

  return(NULL)
}


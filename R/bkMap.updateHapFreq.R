bkMap.updateHapFreq <-
function(bkMap, bkIndex, expression=NULL, freq){
  bkFrame = bkMap$bks[[bkIndex]]

  if( length(freq)!= nrow(bkFrame)) stop(paste("updateBlockBinaFreq:: freq length doesn't match the ", bkIndex, " block config.", sep=""))

  newfreq = freq/(sum(freq))

  if(!is.null(expression)) {
    if( length(expression)!= nrow(bkFrame)) stop(paste("updateBlockBinaFreq:: expression length doesn't match the ", bkIndex, " block config.", sep=""))

    bkFrame[, bkMap$expCol ] = expression
  }

  bkFrame[, bkMap$probCol ] = newfreq

  bkMap$bks[[bkIndex]] = bkFrame

  return(bkMap)
}


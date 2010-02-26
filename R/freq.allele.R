freq.allele <-
function(ct, order = c("majorHM", "heter", "minorHM")){

  if ( sum(is.na(match(order, c("majorHM", "heter", "minorHM")))) != 0) stop ("Wrong order of count.")
  freq = (2*ct[1]+ct[2])/( 2*sum(ct) )
  freqs = c(freq, 1-freq)
  #print(ct)
  if ( freq>1 ) warning("Allele frequency is greater than 1")
  return(freqs)

}


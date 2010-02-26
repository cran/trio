isDigitAtLociIdx <-
function(hapIdx, digit=1, lociCt, lociIdx, intVec=seq.int(from=1, to=2^(lociIdx-1), by=1)){

  ifD = F
  
  if(digit!=1 & digit!=2) stop(paste("Invalide digit imput: ", digit, sep=""))
  
  leftSide = (hapIdx - intVec)/2^lociIdx + 1

  if(ifD) print(leftSide)
  
  roundLeftSide = as.integer(leftSide)

  if(ifD) print(roundLeftSide)

  integerLeft = leftSide[roundLeftSide == leftSide]

  if(ifD) print(integerLeft)

  one = integerLeft >= 1
  oth = integerLeft <= 2^(lociCt-lociIdx)

  fittedCt = one & oth

  if(ifD) print(fittedCt)

  if(length(fittedCt)>1) stop("Error!. Should left with only one or none choice")

  if(length(fittedCt)==1){
    if(digit==2) fittedCt = !fittedCt
    return(fittedCt)
  }else{
    if(digit==2) {
      return(TRUE)
    }else{
      return(FALSE)
    }
  }
}


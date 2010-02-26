qing.cut <-
function(val, cutPt, cutPt.ordered = T, right.include=T){

  ## !!! HARD CODE !!! -1 means the val >/>= the maximum value in the cutPt
  ## return the matched index
  if(!cutPt.ordered) cutPt = order(cutPt)

  cutPt.ct = length(cutPt)

  if(cutPt.ct<1) stop("Zero length cutPt.")
  if(cutPt.ct<1) stop("Zero length val.")
  
  val.ct = length(val)

  re = rep(NA, length=val.ct)
  
  numFalse = re
  cellMatchedIdx = re

  if(right.include){
    for( i in 1:val.ct){
      cVal = val[i]
      falseMat = cVal > cutPt
      numFalse[i] = sum(falseMat)
    }
  }else{
    for( i in 1:val.ct){
      cVal = val[i]
      falseMat = cVal >= cutPt
      numFalse[i] = sum(falseMat)
    }
  }
    
  cellMatchedIdx [ numFalse==cutPt.ct ] = -1
  cellMatchedIdx [ numFalse!=cutPt.ct ] = numFalse[ numFalse!=cutPt.ct ]+1  

  return(cellMatchedIdx)

}


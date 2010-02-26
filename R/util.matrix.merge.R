util.matrix.merge <-
function(ma, ma2){
  maVec = as.vector(ma)
  maVec2 = as.vector(ma2)
  re = c(maVec, maVec2)
  renames = NULL
  if(is.matrix(ma)) {
    maLen = dim(ma)[2]
    renames = colnames(ma)
  }else{
    maLen = 1
    if(is.null(names(ma))){
      renames = ""
    }else{
      renames = names(ma)
    }
  }
  

  if(is.matrix(ma2)) {
    maLen2 = dim(ma2)[2]
    renames = c(renames, colnames(ma2))

  }else{
    maLen2 = 1
    if(is.null(names(ma2))){
      renames = c(renames, "")
    }else{
      renames= c(renames, names(ma2))
    }
  }    
  if(F)print(renames)
  colNum = maLen + maLen2
  if(F)print(colNum)
  re = matrix(re, ncol=colNum)
  colnames(re) = renames
  
  return(re)
}


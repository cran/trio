listExtractor <-
function(nameList, varList, na.replace = NULL){
  ifD = F

  varNames = names(varList)
  varLen = length(varList)
  exLen = length(nameList)

  m = match(nameList, varNames, 0)
  if(min(m)==0 & is.null(na.replace)){
    stop(paste("Variable(s) not found in the varList. NULL not allowed. Names in varList: ", paste(varNames, collapse=", "), sep=""))
  }

  re = NULL
  for(i in 1:length(m)){
    if(m[i]==0){
      re = c(re, na.replace)
    }else{
      re = c(re, varList[m[i]])
    }
  }
  return (re)
}


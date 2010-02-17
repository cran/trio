util.list.ex <-
function(nameList, varList, na.allow = T, na.replace = NULL){
  ifD = F

  varNames = names(varList)
  varLen = length(varList)
  exLen = length(nameList)

  if(varLen!=exLen) warning(paste("Two lists are of different length. nameList of length(", exLen, "); varList of length(", varLen,").", sep="") )
  ## create a map with key = var name, item = var position 
  varMap = data.frame(item = seq.int(from=1, to=varLen, by=1), key = varNames)

  namePos = rep(0, length=exLen)
  
  tryCatchEnv = new.env(parent=baseenv())
  assign("namePos", namePos, env=tryCatchEnv)
  
  ## based on this map, find the right sequence of positions to extract values
  for ( i in 1:exLen){
    tryCatch({namePos[i] = varMap$item[varMap$key == nameList[i]]},
             error = function(e){
               if(!na.allow){
                 stop(paste("var with name=[", nameList[i], "] not found in the varList.", sep=""))
               }else{
				 a = get("namePos", env=tryCatchEnv)
				 a[i]=0
				 assign("namePos", a, env=tryCatchEnv)
                 #namePos[i] <<- 0
               }
             })

  }

  ## assign default value
  varEx = rep(NA, length=exLen)
  if(!is.null(na.replace)){
    varEx = rep(na.replace, length=exLen)
  }
  
  for ( i in 1:exLen ){
    if(namePos[i]!=0){
      varEx[i]= unlist(varList[namePos[i]])
      if(ifD) print(paste("i=(", i, ") namePos=(", namePos[i], ")."))
    }
  }
    
  return (varEx)
}


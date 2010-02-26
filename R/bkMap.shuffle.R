bkMap.shuffle <-
function(bkMap, bkCt=1, keys.new=NULL, exclude=F){
  keys.ori = bkMap$keys
  if(is.null(keys.new)){
    keys.new =  resample(keys.ori, size=bkCt, replace=F)
  }

  # order are index in the keys.ori
  if(!exclude){
    keys.neworder = match(keys.new, keys.ori)
    leftkeys.neworder = (1:length(keys.ori)) [-keys.neworder]

    ## CHANGE. do not think keep the left over blocks is ever be used.
    #if(!del.left){
    #  keys.neworder = c(keys.neworder, leftkeys.neworder)
    #}else{
      keys.neworder = keys.neworder
    #}
  }else{
    # if want to exclude the keys.new
    keys.neworder = match(keys.new, keys.ori)
    keptkeys.neworder = (1:length(keys.ori)) [-keys.neworder]

    #if(!del.left){
    #  keys.neworder = c(keptkeys.neworder, keys.neworder)
    #}else{
      keys.neworder = keptkeys.neworder
    #}
  }

  ##print(keys.neworder)
  newBkMap = bkMap
  newBkMap$bks = bkMap$bks[keys.neworder]
  newBkMap$bkLens = bkMap$bkLens[keys.neworder]
  newBkMap$keys = keys.ori[keys.neworder]
  newBkMap$snpCt = bkMap$snpCt[keys.neworder]
  newBkMap$snpLen = sum(newBkMap$snpCt)
  newBkMap$alleleCode = bkMap$alleleCode

  if(exclude){
    newBkMap$selBkCt = length(keys.ori) - length(keys.new)
  }else{
    newBkMap$selBkCt = length(keys.new)
  }
  return(newBkMap)
  
}


util.str.splitAndOrNot <-
function(str, split="and", check.space=F){
  split.ex = util.str.splitEx(str, split)
  ## for example "and"
  if(length(split.ex)==1){
    ## check whether no "and" is found
    if (nchar(split.ex[1])==nchar(str)) {
      re = NULL
      return(re)
    }
  }
  ## remove space
  split.ex.rm= unlist(lapply(split.ex, util.str.rmSpacePadder, rm=c("b", "e")))
  items=util.array.rmEmptyStr(split.ex.rm)

  ## if check.space, return NULL if any of the item contains space in it
  if(check.space){
    found=F
    i=1
    while(i<=length(items) & !found){
      tmp.idx = util.char.1stIdx(items[i], " ")
      if(tmp.idx>0) found=T
      i = i + 1
    }
    if(found) return(NULL)
  }
  ## create the list for return
  re = list(op=split, items=util.array.rmEmptyStr(split.ex.rm))
  return(re)
 
}


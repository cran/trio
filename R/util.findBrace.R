util.findBrace <-
function(str){

  #left ==0, right=1
  brace=-1
  idx=0
  
  arr = util.str.2CharArray(str, len=nchar(str))
  l.find =util.char.1stIdx(arr, "(", match=T, is.array=T)
  r.find =util.char.1stIdx(arr, ")", match=T, is.array=T)

  if(l.find==0){
    if(r.find!=0)
      return (list(brace=1, idx=r.find))
    else
      return(NULL)
  }
  if(r.find==0){
    if(l.find!=0)
      return (list(brace=0, idx=l.find))
    else
      return(NULL)
  }
  if(l.find<r.find) return (list(brace=0, idx=l.find))
  if(l.find>r.find) return (list(brace=1, idx=r.find))
  
}


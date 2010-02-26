util.str.replace <-
function(str, replaced, new, replace.all=F){

  pos = util.char.1stStrIdx(str, find=replaced)
  if(pos==0){
    # if not found, return the original
    return(str)
  }

  if(  (pos+nchar(replaced))<=nchar(str)){
    str.left = substr(str, pos+nchar(replaced), nchar(str))
  
    if(replace.all){
      # iteratively replace the left-over part
      str.left = util.str.replace(str.left, replaced=replaced, new=new, replace.all=replace.all)
    }
    str.re = paste(substr(str, 1, pos-1), new, str.left, sep="")
  }else{
    str.re = paste(substr(str, 1, pos-1), new, sep="")
  }
  return(str.re)

}


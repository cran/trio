util.str.splitEx <-
function(str, split, case.sen=T){

  if(!case.sen){
    str = tolower(str)
    split = tolower(split)
  }
  split.re = strsplit(str, split=split)
  split.re=split.re[[1]]

  len.s = nchar(split.re)
  re = split.re[len.s!=0]
  return(re)
}


util.str.2CharArray <-
function(str, len=nchar(str)){

  startSeq = 1:len
  re = lapply(startSeq, FUN = function(inStr, startSeq){
    char = substr(inStr, startSeq, startSeq)
    char
  }, inStr = str)
  re = unlist(re)
  return(re)
}


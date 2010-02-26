qstrsplit <-
function(str, delim, re.1st = F){
  re = unlist(strsplit(str, delim))
  if (length(re)==1){
    # possible not find the delim
    if (nchar(re[1])==nchar(str)){
      # not find the delim
      return(NA)
    }else  if(nchar(re[1])==0){
      # str contain delim itself
      if(re.1st ){
        return(re)
      }else{
        return(NULL)
      }
    }else{
      # str ends with delim
      return(re)
    }
  }else{
    if(nchar(re[1])==0){
      # str starts with delim
      if(re.1st){
        return(re)
      }else{
        return(re[-1])      
      }
    }else{
      return(re)
    }
  }
}


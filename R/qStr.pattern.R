qStr.pattern <-
function(objList,  exam.ext = .3){
  l.len = length(objList)

  class.last = NULL
  len.last = NULL
  names.last = NULL
  i = 1
  search=T
  while( i <= min((round(l.len*exam.ext, 0)+1), l.len) & search ){
    #print(i)
    elm = objList[[i]]
    names = names(elm)
    class = unlist(lapply(elm, FUN=class))
    len = unlist(lapply(elm, FUN=length))

    if(!is.null(names.last)){
      ## see if the names matched
      if(length(names.last)==length(names)){
        all.m = names.last==names
        if(sum(all.m)!=length(names)) search=F
      }else{
        search=F  
      }
      if(length(class.last)==length(class)){
        all.m = class.last==class
        if(sum(all.m)!=length(class)) search=F
      }else{
        search=F  
      }
      if(length(len.last)==length(len)){
        all.m = len.last==len
        if(sum(all.m)!=length(len)) search=F  
      }else{
         search=F  
      }      
      
    }
    names.last=names
    names=NULL
    class.last=class
    class=NULL
    len.last=len
    len=NULL    
    i = i+1
  }

  if(search){
    if (is.null(names.last[1])){
      reStr = paste(c("name", "class", "length"),
                  c("",  class.last[1],  len.last[1]), 
                  sep="=", collapse="; ")
    }else{
      reStr = paste(c("name", "class", "length"),
                  c(names.last[1],  class.last[1],  len.last[1]), 
                  sep="=", collapse="; ")
    }
  }else{
    reStr = NULL
  }
  return(reStr)
}


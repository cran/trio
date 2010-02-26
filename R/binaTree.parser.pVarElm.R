binaTree.parser.pVarElm <-
function(str, type="un", check.format=T){
    ifD = F
    fN = "binaTree.parser.pVarElm:"

    str = util.str.rmSpacePadder(str, rm=c("b", "e"))
    if(ifD) print(qp("binaTree.parser.pVarElm::str=", str, "|begin"))
    if(ifD) print(qp("binaTree.parser.pVarElm::type=", type, "|begin"))
    
    idx = NULL
    bool = NA
    mixed = 0
    
    if(type=="un"){
      # see if starting with 'and' or 'or'or 'not'
      idx.and = util.char.1stStrIdx(str, find="and ")
      idx.or = util.char.1stStrIdx(str, find="or ")
      idx.not = util.char.1stStrIdx(str, find="not ")
      if(idx.not==1){ ## start with not
          mixed = 0
          bool = -1
          idx = util.str.splitAndOrNot(str, "not", check.space=check.format)
          if(!is.null(idx)){
            if(check.format){
              if(length(idx$items)==1){
                return(list(idx=idx$items, bool=bool, mixed=mixed))
              }
            }else{
              if (length(idx$items)>2) stop(paste(fN, "do not allow more than two variables in an element."))
  
              return(list(idx=idx$items, bool=bool, mixed=mixed))
            }
          }
      }

      if(idx.and==1){ ## start with and
          mixed = 1
          bool = 0
          idx = util.str.splitAndOrNot(str, "and", check.space=check.format)

          if(is.null(idx)) return(NULL)

          if (length(idx$items)>2) stop(paste(fN, "do not allow more than two variables in an element."))
  
          return(list(idx=idx$items, bool=bool, mixed=mixed))
      }
      if(idx.or==1){ ## start with or
          mixed = 1
          bool = 1
          idx = util.str.splitAndOrNot(str, "or", check.space=check.format)

          if(is.null(idx)) return(NULL)
          if (length(idx$items)>2) stop(paste(fN, "do not allow more than two variables in an element."))
          return(list(idx=idx$items, bool=bool, mixed=mixed))
      }

      list.and = util.str.splitAndOrNot(str, " and ", check.space=check.format)
      list.or  = util.str.splitAndOrNot(str, " or ", check.space=check.format)

      ## operator is in the middle
      #print(list.and)
      #print(list.or)
      
      if (!is.null(list.and)){
        if (length(list.and$items)>2) stop(paste(fN, "do not allow more than two variables in an element."))
        return(list(idx=list.and$items, bool=0, mixed=2))
      }
      
      if (!is.null(list.or)){
        if (length(list.or$items)>2) stop(paste(fN, "do not allow more than two variables in an element."))
        return(list(idx=list.or$items, bool=1, mixed=2))
      }
    }
    
    if(type=="start"){
      list.and = util.str.splitAndOrNot(str, "and ", check.space=check.format)
      list.or  = util.str.splitAndOrNot(str, "or ", check.space=check.format)

      if (!is.null(list.and)){
        if (length(list.and$items)>1) stop(paste(fN, "do not allow more than two variables without a brace."))
        return(list(idx=list.and$items, bool=0, mixed=1))
      }
      if (!is.null(list.or)){
        if (length(list.or$items)>1) stop(paste(fN, "do not allow more than two variables without a brace."))
        return(list(idx=list.or$items, bool=1, mixed=1))  
      }
    }

    if(type=="end"){
      if(str=="and")
        return(list(idx=NULL, bool=0, mixed=1))
      
      if(str=="or")
        return(list(idx=NULL, bool=1, mixed=1))
      
      list.and = util.str.splitAndOrNot(str, "and", check.space=check.format)
      list.or  = util.str.splitAndOrNot(str, "or", check.space=check.format)

      if (!is.null(list.and)){
        if (length(list.and$items)>1) stop(paste(fN, "do not allow more than two variables without a brace."))
        return(list(idx=list.and$items, bool=0, mixed=1))
      }
      
      if (!is.null(list.or)){
        #print(list.or)
        if (length(list.or$items)>1) stop(paste(fN, "do not allow more than two variables without a brace."))
        return(list(idx=list.or$items, bool=1, mixed=1))  
      }
    }    

    
    if(type=="middle"){
      idx.not = util.char.1stStrIdx(str, find="not ")
      if(idx.not==1){ ## start with not
          mixed = 0
          bool = -1
          idx = util.str.splitAndOrNot(str, "not", check.space=check.format)
          if(!is.null(idx)){
            if(check.format){
              if(length(idx$items)==1) return(list(idx=idx$items, bool=bool, mixed=mixed))
            }else{
              if (length(idx$items)>2) stop(paste(fN, "do not allow more than two variables in an element."))
              return(list(idx=idx$items, bool=bool, mixed=mixed))
            }
          }
      }
      list.and = util.str.splitAndOrNot(str, " and ", check.space=check.format)
      list.or  = util.str.splitAndOrNot(str, " or ", check.space=check.format)

      #print(list.and)
      #print(list.or)
      if (!is.null(list.and)){
        if (length(list.and$items)>2) stop(paste(fN, "do not allow more than two variables in an element."))      
        return(list(idx=list.and$items, bool=0, mixed=2))
      }
      if (!is.null(list.or)){
        if (length(list.or$items)>2) stop(paste(fN, "do not allow more than two variables in an element."))      
        return(list(idx=list.or$items, bool=1, mixed=2))
      }
      if(is.null(list.and) & is.null(list.or)){
        # could be just one element
         if(util.char.1stIdx(str, " ")==0){
           return(list(idx=str, bool=0, mixed=0))
         }
      }
                
    }    

    return(NULL) 
}


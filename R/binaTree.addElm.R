binaTree.addElm <-
function(binaTree,  elm.level, elm.var,  elm.bool){
  # objects for variable only elment for memory reason
  # 1)elm.vlist: a list of vector of variable names, one vector for one element
  # 2)elm.vMa: A matrix for variable only element
  #   matrix col 1:id: list id
  #          col 2:level: level
  #          col 3:bool: boolean term: not(-1) vs 1 ## CHANGED 8JAN09 before:[and(0); or(1); not(-1)]
  
  # two objects for state info
  # curElm: whether current element is a variable only(0) or node(1), or missing(NA)
  # curElmId: the index of current element, or missing(0)
  
  # when the element is a "not" element, the level has to be -1
  ifD = F
 
  if (ifD) print("binaTree.addElm")
  if (ifD) print(qp("elm.level=", elm.level))
  if (ifD) print(qp("elm.var=", elm.var))
  if (ifD) print(qp("elm.bool=", elm.bool))
  if(elm.bool==-1) elm.level=elm.level-1
  
  elm.vlist = binaTree$elm.vlist
  elm.vMa = binaTree$elm.vMa
  
  elm.vlist = c(elm.vlist, list(elm.var))
  
  if(binaTree$curElmId!=0){
    
    elm.lastIdx = max(elm.vMa[,1])
    elm.vMa = rbind(elm.vMa, c(elm.lastIdx+1, elm.level, elm.bool))
  }else{
    ## initialize the empty variables
    elm.lastIdx=0
    elm.vMa[1,]= c(elm.lastIdx+1, elm.level, elm.bool)
  }
  binaTree$elm.vlist = elm.vlist
  binaTree$elm.vMa = elm.vMa

  binaTree$curElm = 0;  #variable element
  binaTree$curElmId = elm.lastIdx+1 # id
  if (ifD) print("after addElm, binaTree")
  if (ifD) print(binaTree)
  return(binaTree)
}


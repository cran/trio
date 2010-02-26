binaTree.closeNode <-
function(binaTree,  elm.level){
  ## if the current element is an element with not as boolean,
  ## it will have have row in the node map

  # objects for node element. (contain node itself or variable only element)
  # 1)elm.nlist: a list of vector, each vector for one element
  #    matrix col 1: idx for node or variable element
  #    matrix col 2: is a node(1) or variable element(0)
  # 2)elm.nMa: A matrix for node
  #    matrix col 1:id: list id
  #           col 2:level: levle
  #           col 3:bool: boolean term: and(0); or(1)
  #           col 4:open: open(1)/close(0)
  # two objects for state info
  # curElm: whether current element is a variable only(0) or node(1), or missing(NA)
  # curElmId: the index of current element

  ifD = F
  if (ifD) print("binaTree.closeNode")
  if (ifD) print(qp("elm.level=", elm.level))
  if (ifD) print("binaTree::")
  if (ifD) print(binaTree)
  
  elm.nlist = binaTree$elm.nlist
  elm.nMa = binaTree$elm.nMa
  open.id = NA

  ## print(elm.nMa)
  ## see if any node is open (should be)
  if (nrow(elm.nMa)>=1){
      fil = (elm.nMa[,2]==elm.level) & (elm.nMa[,4]==1)
      if (ifD) print(fil)
      if(sum(fil, na.rm=T)==1){
        # close the old node
        open.id = elm.nMa[,1][fil]
        elm.nMa[fil, 4]=0

        nMa = elm.nlist[[which(fil)]]
        if(ifD) print(nMa)
        nMa = rbind(nMa, c(binaTree$curElmId, binaTree$curElm))
        elm.nlist[which(fil)]=list(nMa)
      }else{
        stop("no node match the requirement")
      }
  }else{
      stop("no open node")
  }
  
  binaTree$elm.nlist = elm.nlist
  binaTree$elm.nMa = elm.nMa
  binaTree$curElmId = open.id#
  binaTree$curElm = 1 #node

  if(ifD) print("after close node, binaTree:")
  if(ifD) print(binaTree)
  
  return(binaTree)
}


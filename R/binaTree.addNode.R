binaTree.addNode <-
function(binaTree,  elm.level, elm.bool){
  ifD = F
  if (ifD) print("binaTree.addNode")
  if (ifD) print(qp("elm.level=", elm.level))
  if (ifD) print(qp("elm.bool=", elm.bool))
   
  # objects for node element. (contain node itself or variable only element)
  # 1)elm.nlist: a list of matrix, each vector for one element
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
  
  ## called when the node is open or adding element
  elm.nlist = binaTree$elm.nlist
  elm.nMa = binaTree$elm.nMa

  ## if elm.nMa has the node open
  if (length(elm.nlist)>=1){
      fil = elm.nMa[,2]==elm.level & elm.nMa[,4]==1
      if(sum(fil, na.rm=T)==1){
        # append to old node
        open.id = elm.nMa[,1][fil]
        
        nMa = elm.nlist[[fil]]
        nMa = rbind(nMa, c(binaTree$curElmId, binaTree$curElm))
        elm.nlist[fil]=list(nMa)
      }else{
        # new node
        elm.lastIdx = max(elm.nMa[,1])
        open.id = elm.lastIdx + 1
        elm.nMa = rbind(elm.nMa, c(open.id, elm.level, elm.bool, 1))

        nMa = matrix(c(binaTree$curElmId, binaTree$curElm), nrow=1, ncol=2)
        elm.nlist=c(elm.nlist, list(nMa))
      }
  }else{
      ## initialize the empty variables
      open.id = 1
      elm.nMa[1,]=c(open.id, elm.level, elm.bool, 1)
      nMa = matrix( c(binaTree$curElmId, binaTree$curElm), nrow=1, ncol=2)
      elm.nlist = c(elm.nlist, list(nMa))
  }
  
  binaTree$elm.nlist = elm.nlist
  binaTree$elm.nMa = elm.nMa
  
  binaTree$curElm = 1;  #node element
  binaTree$curElmId = open.id # id
  if (ifD) print("after addNode, binaTree:")
  if (ifD) print(binaTree)
  return(binaTree)
}


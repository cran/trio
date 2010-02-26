binaTree.shuffle.findNodePar <-
function(tree, id.node){
  fN = "binaTree.shuffle.findNodePar:"
  elm.nMa = tree$elm.nMa

  par.level = elm.nMa[elm.nMa[,1]==id.node,2]-1
  if(par.level<0) {
    #warning(paste(fN, " top level has no parent for id =", id.node))
    return(0)
  }

  par.posiId = elm.nMa[elm.nMa[,2]==par.level,1]

  if(length(par.posiId)<=0 ) stop( paste(fN, " no node belongs to the upper level"))

  i = 1

  posId = NULL
  ## multiple parents may have it as a child
  search = T
  while(i<=length(par.posiId) & search){

    posParId = par.posiId[i]
    nodeMa = tree$elm.nlist[[ posParId ]]
    yesNode = nodeMa[ ,2]==1
    childId = nodeMa[yesNode,1]
    if(sum(is.element(id.node, childId))==1){
      posId = posParId
      search = F
    }
    i = i +1;
  }

  return(posId)
 
}


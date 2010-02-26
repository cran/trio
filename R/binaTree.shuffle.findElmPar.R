binaTree.shuffle.findElmPar <-
function(tree, id.elm){
  fN = "binaTree.shuffle.findElmPar:"
  elm.nMa = tree$elm.nMa

  ## loop through each node in nMa
  search = TRUE
  i = 1
  child.node=NULL
  while ( i <= length(tree$elm.nlist) & search){
    curMa = tree$elm.nlist[[i]]
    filter = curMa[ ,2]==0 & curMa[,1]==id.elm
    if(sum(filter)>=1){
      search = FALSE
      child.node = i
    }
    i = i+1
  }

  if(search){
    stop(paste(fN, " cannot find the node using the current element."))
  }

  return(child.node)

}


binaTree.parser <-
function(str){

  ifD = F

  ## the first character is "(", "not", "and", "or", or variable name
  b.lev=0
  b.open = F
  cur.idx = 0
  cur.str = str

  binaTree = binaTree.constr()

  binaTree = binaTree.parser.proc(str, binaTree, curLevel=0, first.seg=T)

  # check wether the level 0 is closed
  nMa = binaTree$elm.nMa

  if (length(binaTree$elm.nlist)>0){
    cur.level = nMa[, 2]
    if(sum(cur.level==0)>1) stop()
    zero.level = which(cur.level==0)
    if( nMa[zero.level, 4]==1){
      ## close the last node
      binaTree=binaTree.closeNode(binaTree,  elm.level=0)
    }
    
  }
  
  return(binaTree)

}


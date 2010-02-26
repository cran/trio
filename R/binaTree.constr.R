binaTree.constr <-
function(){

  elm.vlist = list()
  elm.vMa = matrix(NA, ncol=3, nrow=1)
  colnames(elm.vMa) = c("id", "level", "bool")

  elm.nlist = list()
  elm.nMa = matrix(NA, ncol=4, nrow=1)
  colnames(elm.nMa) = c("id", "level", "bool", "open")

  binaTree = list(elm.vlist=elm.vlist, elm.vMa=elm.vMa, elm.nlist=elm.nlist, elm.nMa=elm.nMa,
               curElm=NA, curElmId=0)
  return(binaTree)
}


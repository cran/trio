binaTree.1Level.changeSignal <-
function(bina.col, model.signal.cur ){

    bina.col.num = sapply(bina.col, FUN=function(item) {as.numeric(substr(item, 1, nchar(item)-1))})
    bina.col.ct = table(bina.col.num)
    bina.snp.single = paste(dimnames(bina.col.ct)[[1]][bina.col.ct==1], "b", sep="")

    bTree = binaTree.parser(str=model.signal.cur)
    leaves.all = unlist(bTree$elm.vlist)
    leaves.comb = unique(leaves.all[  is.element(unlist(bTree$elm.vlist), bina.snp.single)  ])
    if(length(leaves.comb)>=1){
      for (i in 1:length(leaves.comb)){
        tmp.str = leaves.comb[i]
        ##print(tmp.str)
        model.signal.cur = util.str.replace(str=model.signal.cur,
                           replaced=tmp.str,
                           new=paste(substr(tmp.str, 1, nchar(tmp.str)-1), "a", sep=""),
                           replace.all=T)
        ##print(model.signal.cur)
      }
      
    }
    return(model.signal.cur)
  }


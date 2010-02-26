bkMap.updateByRule <-
function(bkMap, rule, complement=F, is.hapGenoMap = F){

  # shuffle the genotype blocks, but keep the columns
  # so need to shuffle the column names

  if(is.hapGenoMap) stop("Method not implemented for hapGenoMap listing yet!")

  # based on the snp idx, map them with block idx
  if(!is.hapGenoMap){
    causalBkIdx.new = bkMap.findHRCBIdx(bkMap, as.integer(rule$snpIdx), re.keys=F, unique=T, complement=complement)
    total.bkct = length(bkMap$keys)
  
    # construct the beginning/ending index of snp in original map
    idx.be = matrix(NA, nrow=total.bkct, ncol=2)
    idx.be[,2]=cumsum(bkMap$snpCt)
    idx.be[,1]=c(1, cumsum(bkMap$snpCt)+1)[1:total.bkct]
  
    id.be.shuffled = idx.be[ c(causalBkIdx.new,  (1:total.bkct)[!is.element(causalBkIdx.new, 1:total.bkct)]),, drop=F ]
  
    newData.col = apply(id.be.shuffled, 1, FUN= function(item){
      ## use recycling rules
      re = paste("v", rep(item[1]:item[2], each=2), sep="")
      re = paste(re,  c("a", "b"), sep="")
      re})

  }else{
    ## not implemented. Because later we added the argu, complement. 
    ## causalBkIdx.new = hapBkGenoMap.findHRCBIdx(allmap=bkMap, snpIdx=as.integer(rule$snpIdx), re.key=F, unique=T)

    ## snp.idxSeq = hapBkGenoMap.findHRCBSnpIdx(allmap=bkMap, snpIdx=as.integer(rule$snpIdx), re.1digit=T)
    ## newData.col = t( paste( "v", rep(snp.idxSeq, each=2), c("a", "b"), sep=""))

  }

  info = list(causalBkIdx = causalBkIdx.new,
            new.colname = unlist(newData.col))
  return(info)
  

}


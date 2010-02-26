bkMap.findHRCBIdx <-
function(bkMap, snpIdx, re.keys=T, unique=T, complement=F){
  if(!is.integer(snpIdx)) stop()
  # find out the bk bin the snp fall into by comparing the idx with the ending idx of block
  endIdx = cumsum(bkMap$snpCt)
  idx.bk.belong = qing.cut(val=snpIdx, cutPt=endIdx, cutPt.ordered = T, right.include=T)

  if ((max(idx.bk.belong<0)==1) | (max(idx.bk.belong>length(endIdx))==1)) stop("One of the SNP index in sigStr is out of range.")

  if(unique){
    idx.bk.belong = sort(unique(idx.bk.belong)) 
  }
  if(complement){
    total = 1:(length(bkMap$keys))
    idx.bk.belong.unique =  sort(unique(idx.bk.belong)) 
    idx.bk.belong = total[ -idx.bk.belong.unique ]
  }
  
  if(re.keys){
    idx.bk.belong = bkMap$keys[idx.bk.belong]
  }
  return(idx.bk.belong)
}


bkMap.genoFreq <-
function(bkMap, keys){

  # only one key
  idx = match(keys, bkMap$keys)

  genoFreq = NULL

  for( i in idx){
    if( is.na(i) ) stop()
    if(i>1){
      snpBIdx = sum(bkMap$snpCt[ 1:(i-1) ])
    }else{
      snpBIdx = 0
    }
    bks = bkMap$bks[[ i ]]
    genoFreq.byBk = hap.2geno( as.character(bks[, bkMap$expCol]), bks[,bkMap$probCol], snpBIdx)
    genoFreq = rbind(genoFreq, genoFreq.byBk)
  }
  return(genoFreq)

}


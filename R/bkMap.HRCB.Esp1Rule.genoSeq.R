bkMap.HRCB.Esp1Rule.genoSeq <-
function(bkMap, rule, re.probOnly = T){

  ifD = F
  signalRule = rule$slist[[1]]

  HRCBIdx = bkMap.findHRCBIdx(bkMap, as.integer(signalRule$var.idx), re.keys=F, unique=T)

  causalInfo = bkMap.updateByRule(bkMap, rule)
  
  # after process signalRule, 
  keys.new = bkMap$keys[c(causalInfo$causalBkIdx)]
  bkMapS = bkMap.shuffle(bkMap, bkCt=NULL, keys.new=keys.new, exclude=F)

  superHRCBMap = bkMap.superHRCB(bkMap, uniBkIndexes=HRCBIdx, re.probOnly=re.probOnly)

  if(re.probOnly){
    supHapProb = superHRCBMap
    #print(paste("# of prob: ", length( superHRCBMap)))

    return(list(bkMapS=bkMapS, supHapProb=supHapProb))
  }else{
   #print("Re additional stuff")
   #print(superHRCBMap)
   supHapProb = superHRCBMap$prob
   superHapExp = apply(superHRCBMap$linkBkExp, 1, paste, collapse="")

   if(ifD) print( cbind(superHapExp, superHRCBMap$prob, cumsum(superHRCBMap$prob)))
  
   # obtain super-hap and their probabilities

   #print(paste("# of prob: ", length( superHRCBMap$prob)))

   return(list(bkMapS=bkMapS, supHapExp=superHapExp, supHapProb=supHapProb))

  }

}


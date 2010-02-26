bkMap.HRCB.LRCB.split <-
function (bkMap, rule){

  causalInfo = bkMap.updateByRule(bkMap=bkMap, rule=rule)
  col.shuffled=causalInfo$new.colname

  nocausalInfo = bkMap.updateByRule(bkMap=bkMap, rule=rule,  complement=T )
  col.shuffled.nocausal =nocausalInfo$new.colname

  
  # after process signalRule, 
  keys.new = bkMap$keys[c(causalInfo$causalBkIdx)]
  
  # one map for associated blocks
  bkMapS  = bkMap.shuffle(bkMap, bkCt=NULL, keys.new=keys.new, exclude=F)


  bkMapNS = bkMap.shuffle(bkMap, bkCt=NULL, keys.new=keys.new, exclude=T)
        
  return(bkMaps = list(bkMapS=bkMapS,  bkMapNS=bkMapNS, newColName=col.shuffled, noCausalColName=col.shuffled.nocausal))
}


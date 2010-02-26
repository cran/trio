binaTree.patternForm <-
function(binaTree, elmMList, bkIdx.RuleOrder){
  ifD = F
  fN = "binaTree.patternForm:"

  if(ifD) {
    print("binaTree")
    print(binaTree)
    print("elmMList")
    print(elmMList)
    print(bkIdx.RuleOrder)

  }


  # place to hold the matching hap pairs
  mgrid = NULL

  unique.bk = unique(bkIdx.RuleOrder)
  unique.bk = sort(unique.bk)
  
  # total number of HRCB
  bkCt = length(unique.bk)

  nMa = binaTree$elm.nMa
  nCt = nrow(binaTree$elm.nMa)

  ## in binaTree, all the op among node is "OR".
  ## but within a node, it is one element or a boolean term of multiple element
  ## with more than one element.

  mgrid = NULL
  ##process node with "and" op
  filter = nMa[,3]==0
  # find out the node with "and" as operator
  range = (1:nCt)[filter]
  ## need to parse the node with "and" as operator
  if(sum(filter)>=1){
    for ( i in range){
      # for each node 
      node.idx = nMa[i,1]
      nodeMa = binaTree$elm.nlist[[node.idx]]
      if(sum(nodeMa[,2])>=1) stop(paste(fN, "'and' node with node as element."))

      # nodeMa[,1] gives the element's id, which follow the elmBkIdx
      bkIdx.seqInRule = nodeMa[,1]

      bkm.ct = table(bkIdx.seqInRule)
      bkm.bench = NULL
      bkm.bench = as.integer(dimnames(bkm.ct)[1][[1]])

      ## if two elements in one bk need to be meet, need to update elmMList and bkIds
      newMList = NULL
      
      for(j in 1:length(bkm.bench)){
        ## for each bkIdx related to the element in the node
        tmp.ct = bkm.ct[j]
        elm.m.ids =  bkm.bench[j]
        if(tmp.ct>1){ # if more than one element from the same bk
          stop (paste(fN, "not implemented"))
#           elm.m.ids = nodeMa[ bkm.idx==bkm.cur.idx, 1]
#           #obtain all the elmMlist for this set of elm
#           tmp.allMatch = elmMList[elm.m.ids]
#           ## keep the first set of meeting one
#           tmp.finalSet = tmp.allMatch[[1]]
#           final.fit = rep(1, times=nrow(tmp.allMatch))
#           for(ii in 2:tmp.ct){
#             ## do an intersection operation, first assume all are in
#             ## then match two indexes
#             tmp.cur = tmp.allMatch[[ii]]
#             filter = is.element(tmp.finalSet[,1], tmp.cur[,1])
#             ## match the first idx
#             final.fit = final.fit & filter
#             ## match the second idx
#             filter = match(tmp.finalSet[,2], tmp.cur[,2])
#             final.fit = final.fit & filter            
#           } # for(ii in 2:tmp.ct){
#           if(sum(filter)<=0) stop (paste(fN, "no fitted hap pairs for the more than two matching elment in a block."))
# 
#           tmp.ma = tmp.finalSet[filter,, drop=F ]
#           newMList = c(newMList, list(tmp.ma))
        }else{
          newMList = c(newMList, elmMList[elm.m.ids])
        } # if(tmp.ct>1){ # one element for one bk
        
      } # for(j in 1:length(bkm.bench))
      ## form the pattern for add node, add to the match.grid1 and 2
      bkMVec = bkIdx.RuleOrder[ bkm.bench]
      # match the bkm.cur.idx to a position in the whole pattern
      mat.idx = match(bkMVec, unique.bk)
      hapPair.matched=hapPair.match(mat.idx, bkCt, newMList)
      if(ifD) print(hapPair.matched)
      mgrid = rbind(mgrid, hapPair.matched)
    } # for ( i in range){ parse node with "and" as operator
  } # if(sum(filter)>=1){

  ## parse all the element in "or" node

  for (i in 1:nCt){
     curNode = nMa[i,]
     ## only for "or" node
     if(curNode[3]==1){
      node.idx = curNode[1]
      nodeMa = binaTree$elm.nlist[[node.idx]]
      ## only considering the element variables
      if(sum(nodeMa[,2])<nrow(nodeMa)){
        if(ifD) print("find the node linked with one element!!")
        elm.ids = nodeMa[nodeMa[,2]==0,1]
        #if(length(elm.ids)!=1) stop(paste(fN, "number of element in a node is not one"))

        # it is OK to have more than one element
        for (ij in 1:(length(elm.ids))){
          bkm.idx = bkIdx.RuleOrder [ elm.ids[ij] ]
          mat.idx = match(bkm.idx, unique.bk)
          if (ifD) print(bkm.idx)
          if (ifD) print(unique.bk)
          #print(elmMList[elm.ids[ij]])
          hapPair.matched=hapPair.match(mat.idx, bkCt, elmMList[elm.ids[ij]])
          mgrid = rbind(mgrid, hapPair.matched)
        }
        
        if(ifD) print(mgrid)
      }# else it is a all node node      
     } # if(curNode[3]==1){
  } # for i in (1:nCt){
  
  return( mgrid )
}


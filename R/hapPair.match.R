hapPair.match <-
function(listNewOrder, bkMax, pairList){
  ifD=F
  fN = "hapPair.match:"
  if (ifD)   print(paste(fN, " start"))
  if (ifD)   print(listNewOrder)
  if (ifD)   print(bkMax)
  if (ifD)   print(pairList)
    
  bkCt = length(listNewOrder)

  map.row = 2^(bkCt-1)
  if(bkCt ==1) {
    row.ct = nrow(pairList[[1]])
    mgrid1 = matrix(NA, nrow=row.ct, ncol=bkMax)
    mgrid2 = matrix(NA, nrow=row.ct, ncol=bkMax)

    mgrid1[, listNewOrder]= pairList[[1]][,1]
    mgrid2[, listNewOrder]= pairList[[1]][,2]
    #print(cbind(mgrid1=mgrid1, mgrid2=mgrid2))
    return(cbind(mgrid1=mgrid1, mgrid2=mgrid2))
  }else{
    row.ct = cumprod( unlist(lapply(pairList, FUN=nrow)) )[bkCt]
  }

  # map of index for which hap to pull: each column for each index in bkIds,
  idx.map = exhaustHapExp(lociCt=bkCt)$hap
  idx.map = idx.map[1:(nrow(idx.map)/2),,drop=F]

  if(ifD) print(paste("map.row=", map.row, " row.ct=", row.ct))
  
  mgrid1 = matrix(NA, nrow=row.ct*map.row, ncol=bkMax)
  mgrid2 = matrix(NA, nrow=row.ct*map.row, ncol=bkMax)
  
  for(i in 1:map.row){
    if(ifD) print(i)
    # for each row in the map, grap the hap list and do an expand
    grab.hapSet = idx.map[i,]
    ex.list = lapply(1:bkCt, FUN=function(ii, maps, pList){
      ma = pList[[ii]] 
      re = ma[,maps[ii]]
      return(re)
    }, maps=grab.hapSet, pList = pairList)
    ex.grid = qExpandTable(listOfFactor =ex.list, removedRowIdx=NULL, re.row=F)
    if(ifD) print(ex.grid)
    
    row.seq = ((i-1)*row.ct+1) : (i*row.ct) 
    mgrid1[ row.seq,listNewOrder ] = ex.grid

    grab.hapSet = 3-idx.map[i,]
    ex.list = lapply(1:bkCt, FUN=function(ii, maps, pList){
      ma = pList[[ii]] 
      re = ma[,maps[ii]]
      return(re)
    }, maps=grab.hapSet, pList = pairList)
    if(ifD) print(ex.list)
    ex.grid = qExpandTable(listOfFactor =ex.list, removedRowIdx=NULL, re.row=F)
    if(ifD) print(ex.grid)
    
    mgrid2[ row.seq,listNewOrder ] = ex.grid

  }
  #print(cbind(mgrid1=mgrid1, mgrid2=mgrid2))
  return(cbind(mgrid1=mgrid1, mgrid2=mgrid2))

}


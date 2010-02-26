genMatingTBCondOnChild3 <-
function(appVarNames, child, par1, par2, prob1, prob2, logF = NULL, job=1){

  fStr ="[genMatingTBCondOnChild3]:"
  ifD = F
  maxTbl = NULL
  if(ifD) print(par1)
  if(ifD) print(par2)
  tryCatch({
    maxTbl = get(appVarNames$tbl, env=.GlobalEnv)
    maxRow = nrow(maxTbl)
  }, error=function(e){
    errTrace = paste(e, collapse=";", sep="")
    stop(paste("\n", fStr, errTrace, "\nApp-wise Global Variable ", appVarNames$freqMap, " does not exisit."))
  })
  
  tryCatch({
    maxMateTbl = get(appVarNames$mateTbl, env=.GlobalEnv)
  }, error=function(e){
    errTrace = paste(e, collapse=";", sep="")
    stop(paste("\n", fStr, errTrace, "\nApp-wise Global Variable ", appVarNames$freqMap, " does not exisit."))
  })

  childPairNo = nrow(child)

  if(is.null(childPairNo)) {
    child = matrix(child, ncol=2)
    childPairNo = 1
  }

  if(is.null(dim(par1))){
    par1 = matrix(par1, ncol=2)
  }
  if(is.null(dim(par2))){
    par2 = matrix(par2, ncol=2)
  }
  #print("OOOOOLLLDDDDD")
  
  #print(par1)
  #print(par2)
  
  ## get the most probable ones, order the pair by their prob
  par1 = par1[order(prob1,decreasing=T),,drop=F]
  par2 = par2[order(prob2,decreasing=T),,drop=F]

  prob1 = sort(prob1,decreasing=T)
  prob2 = sort(prob2,decreasing=T)
  
  #print(par1)
  #print(par2)
  #print(prob1)
  #print(prob2)
  
  ## would assume that all pairs with different indexes will come at the beginning of the list
  ## and pairs with same indexes will come at the end of the list

  par1short = nrow(par1)<=nrow(par2)

  if(par1short){
    baseNo = nrow(par1)
    othNo = nrow(par2)
    basePair = par1
    othPair = par2
    baseCol = 1:2
    othCol = 3:4
    mateCol = 1
    
  }else{
    baseNo = nrow(par2)
    othNo = nrow(par1)
    basePair = par2
    othPair = par1
    baseCol = 3:4
    othCol = 1:2
    mateCol = 2
  }

  probVec = rep(NA, times=maxRow)
  counter = 0
  ## iterate choice of one parent and the child
  for ( i in 1:baseNo){
    if(ifD) print(paste("i=", i, " out of ", baseNo))
    
    childMeet = apply(child, 1, FUN = function(rowItem, bench){
        re = as.logical(max(is.element(rowItem, bench)))
    }, bench = basePair[i,])

    childSelRowIdx = NULL
    if( sum(childMeet) >=1 ){
        childSelRowIdx = (1:childPairNo)[childMeet]
        for(j in childSelRowIdx){
             ## if two hap in child is the same
             sameHapChild = child[j,1]==child[j,2]
             if(sameHapChild) {
                par2Meet = apply(othPair, 1, FUN = function(rowItem, bench){
                   re = as.logical(max(is.element(rowItem, bench)))
                }, bench = child[j,1])               

             }else{

                childMatch = is.element(child[j,], basePair[i,])
                if(sum(childMatch)==2){
                  othHapReq = child[j,]
                }else if(sum(childMatch)==1){
                  othHapReq = child[j,][!childMatch]
                }else{
                  stop ("programming error")
                }
                if(ifD) print(childMatch)

                par2Meet = apply(othPair, 1, FUN = function(rowItem, bench){
                   re = as.logical(max(is.element(rowItem, bench)))
                }, bench = othHapReq)                  
             }

             if(ifD) print(par2Meet)
             if( sum(par2Meet) >=1 ){
               par2SelRowIdx = (1:othNo)[par2Meet]
               
               if(ifD)  print(paste("par2SelRowIdx=", par2SelRowIdx))
               addedRow = length(par2SelRowIdx)
               if(ifD) print(paste("addedRow=", addedRow))
               if( addedRow != 0){
                  if( (counter + addedRow) > maxRow ) {
                       #rm(maxTbl)
                       #gc()
                       # stop (paste("\n", fStr, "List length (", counter + addedRow , ") exceed the maximum (", maxRow, ") for i =", i, se
                       simpleResult = matTBCondOnChild3.finalProc(maxTbl=maxTbl, counter=counter, maxMateTbl=maxMateTbl,
                               probVec=probVec, prob1=prob1, prob2=prob2, job=job, logF=logF  )
                       return(simpleResult)

                  }
                  tmpPrb = rep(.25, times=addedRow)
                  othSameHap = NULL
                  othSameHap = othPair[par2SelRowIdx,1]==othPair[par2SelRowIdx,2]
  
                  if(ifD)  print(sum(othSameHap))
                  if(sum(othSameHap)>0) {
                    # if other parent are homo, double =.5
                    tmpPrb[othSameHap]=tmpPrb[othSameHap]*2
                  }
                  # if base parent are homo, double for the second time =1
                  if(basePair[i,1]==basePair[i,2]) tmpPrb = tmpPrb*2

                  ## depend on the child is heto, same pairs for parents and child has different prob than different pairs for parents
                  ## parent must be same heto, homo kid will have .25, heter kid will have .5 
                  if( (sum(othSameHap)<addedRow) && (basePair[i,1]!=basePair[i,2]) && ( child[j,1]!=child[j,2])  ){
                    ## 1st rule, there are some hapPairs not the same for the other parent
                    ## 2nd rule, hapPairs for base parent are not the same
                    ## 3nd rule, hapPairs for child are not the same 
                    tmpMaRow = (1:addedRow)[!othSameHap]
                    othPairMa = othPair[par2SelRowIdx[tmpMaRow],,drop=F]

                    ## compare the hapPairs for the other parents with the hapPairs for the base parents.
                    rowTmpMeet = apply(othPairMa, 1, FUN = function(rowItem, bench){
                      re = as.logical(sum(is.element(rowItem, bench))==2)
                    }, bench = basePair[i,])
                    if(sum(rowTmpMeet)>0) tmpPrb[tmpMaRow[rowTmpMeet]]=.5
                  }
                  
                  probVec[(counter+1):(counter+addedRow)]= tmpPrb
                  
                  ##print(matrix(.Internal(rep(par1[i,], addedRow)), ncol=2, byrow=T))
                  ##print( matrix(par2[par2Meet, ], ncol=2))
                  maxTbl[(counter+1):(counter+addedRow), baseCol] = matrix(rep(basePair[i,], addedRow), ncol=2, byrow=T)
                  maxTbl[(counter+1):(counter+addedRow), othCol] = othPair[par2SelRowIdx, , drop=F]
                  ## newly added
                  maxTbl[(counter+1):(counter+addedRow), 5:6] = matrix(rep(child[j,], addedRow), ncol=2, byrow=T)

                  if(ifD) print(maxTbl[(counter+1):(counter+addedRow),])
                  maxMateTbl[(counter+1):(counter+addedRow), mateCol] = rep(i, addedRow)
                  maxMateTbl[(counter+1):(counter+addedRow), 3-mateCol] = par2SelRowIdx
         
                  counter = counter + addedRow
                  if(ifD)print(paste("counter=", counter))
               }
           } ## if( sum(par2Meet) >=1 ){

        } ## for(j in 1:childSelRowIdx){
    } ## if( sum(childMeet) >=1 ){
  }
  simpleResult = matTBCondOnChild3.finalProc(maxTbl=maxTbl, counter=counter, maxMateTbl=maxMateTbl,
                               probVec=probVec, prob1=prob1, prob2=prob2, job=job, logF=logF  )
  return(simpleResult)
}


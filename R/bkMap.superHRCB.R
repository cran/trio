bkMap.superHRCB <-
function(bkMap, uniBkIndexes=c(1:length(bkMap$keys)), re.probOnly = F){

  ifD = F
  ## find every block expression for every block
  linked = lapply(uniBkIndexes, FUN=function(index, bkMap){
                                            bkMap$bks[[index]]}, bkMap=bkMap)

  hapLen = bkMap$snpCt
  
  len = length(uniBkIndexes)
  if(ifD) print(paste("block len=", len, sep=""))
  
  gridProb = lapply(linked, FUN=function(bk, probCol){
                                 re = bk[, probCol]}, probCol=bkMap$probCol)


    ## only calculate the probability
     gridProb = lapply(linked, FUN=function(bk, probCol){
                                   re = bk[, probCol]}, probCol=bkMap$probCol)
    if(length(uniBkIndexes)==1){
       prob = as.matrix(gridProb[[1]])
    }else{
      prob = qExpandTable(listOfFactor = gridProb)
      prob = apply(prob, 1, prod)
    }
    prob = prob/sum(prob)


  if(re.probOnly){
    return(prob)
  }
  
  gridBase = lapply(linked, FUN=function(bk, expCol){
                                 re = as.character(bk[, expCol])}, expCol=bkMap$expCol)
  gridBaseIndex = lapply(bkMap$bkLens[uniBkIndexes], FUN=function(bkLen){
                                 re = 1:bkLen})
  
  if(length(uniBkIndexes)==1){
     mas = as.matrix(gridBase[[1]])
     mashIn = as.matrix(gridBaseIndex[[1]])
  }else{
    mas = qExpandTable(listOfFactor =  gridBase )
    mashIn = qExpandTable(listOfFactor =  gridBaseIndex )

  }

  if(ifD) {
    print(paste("mas index dim=", ""))
    print(str(mas))
    print(mas[1:2,])
    print(mashIn[1:2,])
    print(prob[1:2])
  }
  #updateBkMap = specifyBks(oriData, uniBkIndexes, bkMap, keyCol, hapLenCol)
  #print(mas)

  re = list(linkBkExp=mas, linkBkIn = mashIn, prob=prob)

  ## print(mashIn)
  return(re)  
}


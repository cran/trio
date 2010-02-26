qExpandTable <-
function(listOfFactor =  list( 1:3, 10:13), removedRowIdx=NULL, re.row=F ){
   tblDim = length(listOfFactor)
   tblDimSeq = unlist(lapply(listOfFactor, FUN=length))

   recyMa = matrix(listOfFactor[[1]], ncol=1, nrow=tblDimSeq[1])
   recyRow = tblDimSeq[1]
   ## grow the index list  
   for ( i in 2:tblDim ){
     growing = util.matrix.clone(recyMa, tblDimSeq[i])
     ## addecCol =  rep(listOfFactor[[i]], each=recyRow)
     recyMa = cbind(growing,  rep(listOfFactor[[i]], each=recyRow))
     recyRow = recyRow * tblDimSeq[i]  
   }

   rowIdx = NULL
   if ( (!is.null(removedRowIdx)) | re.row){
     lengthRev = c(0, cumprod(tblDimSeq)) [ tblDim:1 ]
     tmp = ( removedRowIdx[tblDim:1] - 1)* lengthRev
     rowIdx = sum(tmp)+removedRowIdx[1]
     ## print(lengthRev)
     ## print(tmp)
   }

   ## print(rowIdx)
   ## print(recyMa)
   if(re.row){
     return(list(ma=recyMa, rowNum = rowIdx))
   }else{
     if(!is.null(removedRowIdx)){
       recyMa = recyMa[-rowIdx, , drop=F]
     }
     return(recyMa)
   }
}


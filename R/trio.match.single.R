trio.match.single <-
function(caseRow, controlDf, pedCol, memCol, dadCol, momCol){
   dadRow = which(controlDf[,pedCol]==unlist(caseRow[pedCol]) & controlDf[,memCol]==unlist(caseRow[dadCol]))
   momRow = which(controlDf[,pedCol]==unlist(caseRow[pedCol]) & controlDf[,memCol]==unlist(caseRow[momCol]))

   if( length(dadRow)!=1 ){
     warning(paste("\n The individual is affected, but we have no data on the father. Case information::\n",
                   paste(wrComTbl(caseRow[1:7], colName=colnames(caseRow)[1:7]), collapse="\n"), sep=""))
     return(NULL)
   }
   if( length(momRow)!=1 ){
     warning(paste("\n The individual is affected, but we have no data on the mother. Case information::\n",
                   paste(wrComTbl(caseRow[1:7], colName=colnames(caseRow)[1:7]), collapse="\n"), sep=""))
     return(NULL)
   }
  
   re = rbind(controlDf[dadRow,], controlDf[momRow,], caseRow)
   return(re)
 }


util.matrix.catm <-
function(inMatrix, colVec, discIn=NULL, discVec=NULL, delimVec, digitVec, missingVec=NULL){

     anyMatrix = inMatrix[,colVec]
     mRow = dim(anyMatrix)[1]
     mCol = dim(anyMatrix)[2]
     
     re = matrix(ncol =2, nrow = mRow)

     if(!is.null(discIn)){
       re[,1]= inMatrix[,discIn]
     }else{
       re[,1]=discVec
     }
     for(row in 1:mRow) {
       curRow = NULL
       for(col in 1:mCol){
         num = anyMatrix[row, col]
         if(is.na(num)){
           if(is.null(missingVec)){
             num = "NA"
           }else{
             num = missingVec[col]
           } ## if(is.null(missingVec)){           
         }else{
           ##num = round(num, digitVec[col])
           if(digitVec[col]==0){
             num=format(num, nsmall=0)
           }else{
             num = format(num, digits=digitVec[col], nsmall=digitVec[col])
           }
         } ##if(is.na(num)){
         
         if(col == 1){
           curRow = paste(delimVec[(0+col)], num, sep="")
         }else{
           if(col == mCol){
             curRow = paste(curRow, delimVec[col], num, delimVec[1+col], sep="")
           }else{
             curRow = paste(curRow, delimVec[col], num, sep="")
           } ##if(col = mCol){
         } ## if(col = 1){
         
       } ##for(col in 1:mCol){
       re[row,2] = curRow
     } ##for(row in 1:mRow) {
     return(re)
}


util.matrix.csvText <-
function(numericMa, fileName=NULL, colNames=NULL, rowNames=NULL, digit = NULL,  keepZero=F){
     mRow = dim(numericMa)[1]
     mCol = dim(numericMa)[2]
     i = 0

     lineComment = ""
     
     comm = vector( )

            if(!is.null(colNames)){
              if(!is.null(rowNames))     lineComment = " \t"
              for( col in 1:mCol ){
                lineComment = paste(lineComment, colNames[col], sep="\t")  
              }
              i = i+1
              comm[i] = lineComment
            }     

     lineComment = ""
     
     for( row in 1:mRow ){
       i = i+1
          for( col in 1:mCol ){
          
                if(col == 1){
     
                  if(!is.null(rowNames)){
                        lineComment = paste(lineComment, rowNames[row], sep="\t")  
                      }
                  }
                  tmp = numericMa[row,col]
                  if(is.numeric(tmp) & (!is.null(digit))){
                    if(!keepZero){
                      if(digit==0){
                          num = format(tmp, nsmall=0)
                      }else{
                          num = format(tmp, digits=digit, nsmall=digit)
                      }
                      lineComment = paste(lineComment, num, sep="\t")
                    }else{
                      lineComment = paste(lineComment, round(tmp, digit), sep="\t")
                    }
                  }else{
                    lineComment = paste(lineComment, tmp, sep="\t") 
                  }
     
            }  ## for( col in 1:mCol ){
          comm[i] = lineComment
          lineComment = ""
     }  ## for( row in 1:mRow ){

     if(!is.null(fileName)) write(comm, file = fileName,  append = TRUE)
     return(comm)
}


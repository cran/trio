snpPREFileMatchTrio <-
function(txtF, sep=" ", header = FALSE, pedCol=1, memCol=2, affectCol=6, dadCol=3, momCol=4, txt.affect=2, logF = NULL){

      if(is.character(txtF)){
        ## by default the text file has no header 
        df = read.table(file=txtF, header = header, sep=sep)
      }else{
        df = txtF
      }
      filter = df[, affectCol]==txt.affect 
      caseDf = df[filter,]

      ## it is ok for the parent to be affected
      ##controlDf = df[!filter,]
      controlDf = df

      ## don't sort the case, keep the original order
      #caseDf = caseDf[order(caseDf[,pedCol], caseDf[,memCol]),]
      #controlDf = controlDf[order(controlDf[,pedCol], controlDf[,memCol]),]

      numCase = nrow(caseDf)

      trioDf = NULL
      for( i in 1:numCase){
        caseRow = caseDf[i,]
        
        tryCatch({
             trioDf = rbind(trioDf, trio.match.single(caseRow, controlDf, pedCol, memCol, dadCol, momCol))
           },
             warning = function(warn){
                  if(!is.null(logF)){
                    #print(warn)
                    logWarn(logF, as.character(warn))
                  }else{
                    print("Throw warnings by case")
                    print(warn)
                }
        }) ## tryCatch({
      }
      return(trioDf)
}


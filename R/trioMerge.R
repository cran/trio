trioMerge <-
function(trioData1, trioData2, colName1, colName2,snpCoding=0:3){

  ifD = F

  if( sum(snpCoding == (0:3))!=4)
    stop ("Function use 1-digit coding for genotypes as the following, integer 1 to 3 to
                    represent the three genotypes: less common homozygous, common homozygous, and heterzygous.
                     And zero for missing value")

  trioData = cbind(trioData1, trioData2)

  ## note this is 1-digit coding.
  colName1.first = colName1[ seq.int(from=1, to=length(colName1), by=2) ]
  colName2.first = colName2[ seq.int(from=1, to=length(colName2), by=2) ]

  
  trioCol = c(colName1.first, colName2.first)

  id.all = 1:(length(trioCol))
  
  ori.col = paste("v", id.all, "a", sep="")

  shuffle.order = match(ori.col, trioCol)

  #print( cbind(trioCol, ori.col, shuffle.order))

  trioData.shuffled = trioData[, shuffle.order]

  if(ifD) print(str(trioData.shuffled))
  
  return(trioData.shuffled)

}


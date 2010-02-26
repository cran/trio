setTrioMissingSNP <-
function(trioDf, cord, snp1digit=F, missingDigit = 0){
  ## cord has the beginning row number and col numbers for the trio with missing data

  cord = matrix(cord, ncol=4, byrow=F)
  rowCt = nrow(cord)

  dataNew = trioDf
  if(!snp1digit){
    for( i in 1:rowCt){
      x.y = cord[i,1:2]
      dataNew[  x.y[1]: (x.y[1]+2),  x.y[2]: (x.y[2]+1) ] = missingDigit
      #print(paste("i=", i, " x.y=[", paste(x.y, collapse=";", sep=""), "]", sep=""))
      #print(paste(    x.y[1]: (x.y[1]+2), collapse=";", sep=""))
      #print(paste(    x.y[2]: (x.y[2]+1), collapse=";", sep=""))
    }
  }else{
    for( i in 1:rowCt){
      x.y = cord[i,3:4]
      #dataNew[  ((x.y[1]-1)*3+1): (x.y[1]*3),  x.y[2]+2 ] = 0
      #print(paste("i=", i, " x.y=[", paste(x.y, collapse=";", sep=""), "]", sep=""))
      #print(paste( ((x.y[1]-1)*3+1): (x.y[1]*3) , collapse=";", sep=""))
      #print(paste( x.y[2]+2, collapse=";", sep=""))

      dataNew[  ((x.y[1]-1)*3+1): (x.y[1]*3),  x.y[2]+2 ] = missingDigit
    }
  }
  return(dataNew)
}


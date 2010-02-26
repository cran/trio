util.tbl.exFactorInfo <-
function(tbl){
  rnum = dim(tbl)[1]
  cnum = dim(tbl)[2]

  headList = dimnames(tbl)[[2]]
  numList  = dimnames(tbl)[[1]]
  
  num = NULL
  head = NULL
  ct = NULL
  for ( i in 1:rnum){
    for (j in 1:cnum){
      cur = tbl[i,j]

      num = c(num, numList[i])
      head = c(head, headList[j])
      ct = c(ct, cur)
        
    }
  }
  re = data.frame( row.level=num, col.level=head, ct=as.numeric(ct))
  return(re)
}


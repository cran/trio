util.matrix.2list <-
function(ma, byRow=T){

  if(byRow){
    colNum = dim(ma)[2]
    rowNum = dim(ma)[1]
    ma = t(ma)
  }else{
    colNum = dim(ma)[1]
    rowNum = dim(ma)[2]
  }

  maItem = unlist(as.list(ma))
  re = NULL
  for(row in 1:rowNum){
    cur = maItem[seq.int(from=(row-1)*colNum+1, to=row*colNum, by=1)]
    re = c(re, list(cur))
  }
  return(re)
}


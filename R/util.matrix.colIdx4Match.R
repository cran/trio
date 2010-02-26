util.matrix.colIdx4Match <-
function(ma, val){

  rowCt = nrow(ma)

  #print("util.matrix.colIdx4Match")
  #print(dim(ma))

  matchIdx = qing.mulMatch(val, ma)
  ##matchIdx = which (ma == val)
  
  if(length(matchIdx)>0 & matchIdx[1]!=0){
    
    re = qing.cut(matchIdx,
                  cutPt = seq.int(from=rowCt, to=rowCt*ncol(ma), by=rowCt), 
                  cutPt.ordered = T, right.include=T)

    ##  remove because it cause memeory problems
    ##  as.numeric(cut(matchIdx, breaks = c(0, seq(rowCt, rowCt*ncol(ma), by=rowCt)), include.lowest=T, right=T))
    
  }else{
    
    re =  NULL
  }

  return(re)
}


util.vec.2maIdx <-
function(rowCt, idx){
  ## assuming putted in  column by column
  rowIdx = idx%%rowCt
  colIdx = (idx-rowIdx)/rowCt+1
  if(rowIdx==0) 
    rowIdx=rowCt
    
  return(c(rowIdx, colIdx))
}


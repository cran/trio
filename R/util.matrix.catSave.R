util.matrix.catSave <-
function(data, cols, sep ="-"){
  colNum = dim(data)[2]
  keys = util.matrix.cat(data, cols, sep)
  data[,colNum+1]=keys
  return(data)
}


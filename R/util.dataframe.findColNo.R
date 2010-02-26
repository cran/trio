util.dataframe.findColNo <-
function(data, colNameSearched){
  colNms = colnames(data)
  matchIdx = match(colNameSearched, colNms)
  return(matchIdx)
}


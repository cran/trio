util.str.tokenPicker <-
function(str, delim, indexPicked){
  tokens = unlist(strsplit(str, delim))
  return(tokens[indexPicked])
}


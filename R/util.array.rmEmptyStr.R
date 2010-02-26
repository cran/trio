util.array.rmEmptyStr <-
function(vec){
  lens= nchar(vec)
  re = vec[lens>=1]
  return(re)
}


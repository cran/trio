util.vec.replace <-
function(vec, orignal, replaceBy){
     vec.idx = match(vec, orignal)
     re = replaceBy[vec.idx]
     return(re)
}


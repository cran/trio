util.vec.exByKey <-
function(vecKey, vec, keysOrder){
    filter =  match(keysOrder, vecKey)

    ex = vec[filter]

    return(ex)

}


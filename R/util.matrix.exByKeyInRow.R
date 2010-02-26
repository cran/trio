util.matrix.exByKeyInRow <-
function(ma, keyCol, keys){

    maNew = ma[ is.element( ma[,keyCol], keys), ]

    return(maNew)
}


util.list.rmByKeyVal <-
function(dataList, keys, keyRemoved){
    seqOrder = is.element(keys, keyRemoved)
    return(dataList[!seqOrder])
}


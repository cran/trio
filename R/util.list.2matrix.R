util.list.2matrix <-
function(list, byRow = T, add.dimnames=F) 
{
    colCt = length(list)
    rowCt = length(list[[1]])
    re = matrix(unlist(list), ncol = colCt, nrow = rowCt, byrow = F)
    if(add.dimnames){
      header = dimnames( list[[1]] )[[1]]
      if(is.null(header)) header = paste("v", 1:colCt, sep="")
      rownames(re)=header
    }
    if (byRow) 
        re = t(re)
    return(re)
}


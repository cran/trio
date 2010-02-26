util.listMatrix.2matrix <-
function(list, rbind=T){

  rowCt = length(list)

  re = NULL
  if(rbind){
      for ( i in 1:rowCt ){
         re = rbind(re, list[[i]])
      }
  }else{
      for ( i in 1:rowCt ){
        re = cbind(re, list[[i]])
      }
  }

  return(re)
}


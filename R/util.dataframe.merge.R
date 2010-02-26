util.dataframe.merge <-
function(list, vertical=T){

  len = length(list)
  data = NULL
  for( i in 1:len){
    cur = list[[i]]
    if(vertical){
      data = rbind(data, cur)
    }else{
      data = cbind(data, cur)
    }
  }
  return(data)
}


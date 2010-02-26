util.matrix.cat <-
function(data, cols, sep=""){
  len = length(cols)
  for(i in 1:len){
    if(i==1) {
      re = data[,cols[i]]
    }else{
      re = paste(re, data[,cols[i]], sep=sep)
    }
  }
  return(re)
}


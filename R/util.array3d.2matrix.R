util.array3d.2matrix <-
function(arr, dimLevel=dimnames(arr)[[3]], showDim3 = F, re.num=T, appAsRow = F){
  re = NULL   
  for( i in 1:length(dimLevel)){
    ma = arr[,,i]
    if(showDim3){
      if(re.num){
        ma = cbind(ma, rep(as.numeric(dimLevel[i]), nrow(ma)))
      }else{
        mat = data.frame(ma)
        dimnames(mat)=dimnames(ma)
        ma = cbind(mat, rep(dimLevel[i], nrow(ma)))
      }
    }
    if(appAsRow){
      re = rbind(re, ma)
    }else{
      re = cbind(re, ma)
    }
  }
  return(re)
}


qing.mulMatch <-
function(val, bench){
  
  matched = bench==val
  bench.seq = 1:length(bench)
  re = bench.seq[matched]
  if(length(re)>=1){
    return(re)
  }else{
    return(0)
  }
}


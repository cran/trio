filterHaps2SHaps.check <-
function(hapForBk, hapCts){
  ifD = F
  fN = "filterHaps2SHaps.check::"
  hapList = lapply(hapCts, FUN=function(i) {1:i})
  aa = qExpandTable(listOfFactor =hapList, removedRowIdx=NULL, re.row=F)

  filter = !is.na(hapForBk)
  comp = rep(TRUE, times=nrow(aa))

  
  for( i in 1:length(hapForBk)){
    if(!is.na(hapForBk[i])){
      ff = aa[,i]==hapForBk[i]
      comp =  comp & ff
    }
  }
  row.idx = (1:nrow(aa))[comp]
  if(ifD) print(row.idx)
  
  my.idx = filterHaps2SHaps(hapForBk, hapCts)

  dis = sum(is.na(match(my.idx, row.idx)))
  len.m = length(my.idx)-length(row.idx)

  if (sum(dis, len.m)==0) return (NULL)
  print("compare:my.idx with row.idx")
  print(my.idx)
  print(row.idx)
  stop(paste(fN, "ERROR: not matched index."))
  return(NULL)
}


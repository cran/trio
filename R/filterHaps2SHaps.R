filterHaps2SHaps <-
function(hapForBk, hapCts){
  ifD = F
  fN = "filterHaps2SHaps:"
  if(ifD) {
    print(paste(fN, "begin"))
    print(hapCts)
  }
  bkCt = length(hapCts)

  if(bkCt==1) {
    if (is.na(hapForBk)) {
      return(1:hapCts)
    }else{
      return (hapForBk)
    }
  }

  segLen = cumprod(hapCts)
  
  filter = is.na(hapForBk)

  cut.idx = (1:bkCt)[filter]
  endWithNA = TRUE

  # if the whole pattern doesn't end with NA
  if(!filter[bkCt]){
    endWithNA = FALSE
    cut.idx.end = c(cut.idx, bkCt)
    cut.idx.sta = c(1, cut.idx+1)
  }else{
    # if the whole pattern ends with NA
    cut.idx.end = cut.idx
    cut.idx.sta = c(1, cut.idx[-length(cut.idx)]+1) 
  }

  if(ifD) print(cbind(cut.idx.sta, cut.idx.end))
  segCt = length(cut.idx.sta)
  j = 1
  grp.idx = NULL
  
  while( j <=segCt){
    if(ifD) print(paste("j=", j))
    seg = hapForBk[cut.idx.sta[j]:cut.idx.end[j]]
    if(ifD) print(seg)
    offset.q = 0
    segSeq = 0
    ## calculate the offset
    if( (cut.idx.end[j]-cut.idx.sta[j])!=0){
      ## if the seg has other than NA, we have offset
      if(j==1 & j==segCt & !endWithNA){
        hapCc = hapCts[ (cut.idx.sta[j]):(cut.idx.end[j]) ]
        nums = seg
        offset.q = calHapIdx2SHap(nums, hapCts=hapCc)
        return(offset.q)
      }
      if(j==1 & j==segCt & endWithNA){
        hapCc = hapCts[ (cut.idx.sta[j]):(cut.idx.end[j]-1) ]
        nums = seg[ - length(seg)]
        offset.q = calHapIdx2SHap(nums, hapCts=hapCc)
        n.prec = segLen[cut.idx.end[j]-1]
        segSeq = seq.int(from=0, to=segLen[cut.idx.end[j]]-1, by = n.prec)
        return(segSeq+offset.q)
      }      
      if(j==1){
        ## if it is the first segment and with length >= 2, must have an offset
        hapCc = hapCts[ (cut.idx.sta[j]):(cut.idx.end[j]-1) ]
        nums = seg[ - length(seg)]
        offset.q = calHapIdx2SHap(nums, hapCts=hapCc)
      }
      if(j>1 & j< segCt){
        ## if it is the middle segment and with length >= 2, must have an offset
        hapCc = hapCts[(cut.idx.sta[j]):(cut.idx.end[j]-1) ]
        nums = c(1, seg[ - length(seg)])
        hapCc = c(segLen[(cut.idx.sta[j]-1) ], hapCc)
        ## the out of box offset added one already, but we want to add another offset on it 
        offset.q = calHapIdx2SHap(nums, hapCts=hapCc)-1
      }
      if(j==segCt){
        if(endWithNA){
          ## if it is the last segment and with length >= 2 and have a sequence cutoff
          hapCc = hapCts[(cut.idx.sta[j]):(cut.idx.end[j]-1) ]
          nums = c(1, seg[ - length(seg)])
          hapCc = c(segLen[(cut.idx.sta[j]-1) ], hapCc)
          ## the out of box offset added one already, but we want to add another offset on it 
          offset.q = calHapIdx2SHap(nums, hapCts=hapCc)-1
        }else{
          ## if it is the last segment and with length >= 2 and have no sequence cutoff
          hapCc = hapCts[(cut.idx.sta[j]):(cut.idx.end[j]) ]
          nums = c(1, seg)
          hapCc = c(segLen[(cut.idx.sta[j]-1) ], hapCc)
          ## the out of box offset added one already, but we want to add another offset on it 
          offset.q = calHapIdx2SHap(nums, hapCts=hapCc)-1
        } ## if(endWithNA){
      } ## if(j==segCt){
      if(ifD) print(offset.q)

      if(j!=segCt | endWithNA){
        ## it has a sequence as well
        if(cut.idx.end[j]==1){
          n.prec = 1
          
        }else{
          n.prec = segLen[cut.idx.end[j]-1]
        }        
        segSeq = seq.int(from=0, to=segLen[cut.idx.end[j]]-1, by = n.prec)
      }
        
    }else{ ###  if( (cut.idx.end[j]-cut.idx.sta[j])!=0){
      ## if only one element in the segment

      if (j==segCt & !endWithNA){
          ## if it is the last segment and with length >= 2 and have no sequence cutoff
          hapCc = hapCts[(cut.idx.sta[j]):(cut.idx.end[j]) ]
          nums = c(1, seg)
          hapCc = c(segLen[(cut.idx.sta[j]-1) ], hapCc)
          ## the out of box offset added one already, but we want to add another offset on it 
          offset.q = calHapIdx2SHap(nums, hapCts=hapCc)-1
          segSeq = 0
      }else{
        ## it has a sequence only
        if(cut.idx.end[j]==1){
          n.prec = 1
          ## for the first element to be a NA
          offset.q = 1
        }else{
          n.prec = segLen[cut.idx.end[j]-1]
        }        
        segSeq = seq.int(from=0, to=segLen[cut.idx.end[j]]-1, by = n.prec)
      } ## if (j==segCt & !endWithNA){

    } ###  if( (cut.idx.end[j]-cut.idx.sta[j])!=0){

    if(ifD) print("segSeq")
    if(ifD) print(segSeq)
    if(ifD) print(offset.q)

    if(j==1){
      ## first seg
      grp.idx = segSeq + offset.q
    }else{
      oth.idx = segSeq + offset.q
      all.idx = qExpandTable(listOfFactor =list(grp.idx, oth.idx), removedRowIdx=NULL, re.row=F)
      if(ifD) print(all.idx)
      grp.idx = all.idx[,1]+all.idx[,2]
    }
    if(ifD) print(grp.idx)
    j = j+1
  }
  return(grp.idx)

}


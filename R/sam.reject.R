sam.reject <-
function(allIdx=NULL, idxProb=NULL, rejectIdx, size=1){
	
	if (is.null(allIdx)){
		## if allIdx==NULL, then allIdx = 1:length(idxProb)
		length = length
		allIdx = 1:length
	}else{
		## if idxProb==NUL, then prob are the same for all idx,
		length = length(allIdx)
		if(is.null(idxProb)) idxProb = rep(1/length, length)
	}
	
	if(length<=0) stop(paste("Invalid input value for allList with length=[", length, "]", sep=""))
	
	meet = match(rejectIdx, allIdx)
	idxProb[meet]=0
	idx = resample(allIdx, size=size, prob=idxProb, replace=T)
	
	return(idx)
}


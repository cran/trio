sampleIdxOutsideList <-
function(allList, listVec, maxIt = 1000){
	
	length = length(allList)
	if(length<=0) stop(paste("Invalid input value for allList with length=[", length, "]", sep=""))
	
	keepSearch = T
	idx = NULL
	count = 1 
	while(keepSearch & count<= maxIt){
		idx = sample(1:length, size=1)
		if(!is.element(allList[idx], listVec)){
			keepSearch = F
			return(allList[idx])
		}   
		count = count+1
	}
	
	stop(paste("\nInefficient sampling. Iteration count exceeds maximum limit=[", maxIt, "].",
					"\nallList = [", paste(allList, collapse=";", sep=""),
					"] and listVect = [", paste(listVec, collapse=";", sep=""), "].",
					sep=""))
	return(NULL)
}


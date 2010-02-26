sampleDipSemiAugMap <-
function(semiMapFrame, resiProbCol, augIdxCol, probCol, snpLen){
	## SemiAugMap only have the major haplotypes' expression and prob,
	## only a number for prob for each of the other haplotypes 
	## SemiAugMap also has a column showing the index of the haplotype in the fixed AugMap
	ifD = F
	
	leftOverProb = semiMapFrame[1, resiProbCol]
	mapProb = semiMapFrame[ , probCol]
	mappedAugIdx = semiMapFrame[ , augIdxCol]
	
	commonProbIdx = length(mapProb) 
	
	## need to take account the situation where all the hap is presented.
	if(commonProbIdx == 2^snpLen){
		chooseDip1 = sample(1:commonProbIdx, size=1, prob = mapProb)
		chooseDip2 = sample(1:commonProbIdx, size=1, prob = mapProb)
		
		chooseDip1 = mappedAugIdx[chooseDip1]
		chooseDip2 = mappedAugIdx[chooseDip2]
	}else{
		commonProbIdx = length(mapProb) + 1
		chooseDip1 = sample(1:commonProbIdx, size=1, prob = c(mapProb, leftOverProb))
		chooseDip2 = sample(1:commonProbIdx, size=1, prob = c(mapProb, leftOverProb))
		
		allIdx = NULL
		
		if(chooseDip1 == commonProbIdx){
			allIdx = 1 :(2^snpLen)
			chooseDip1 = sampleIdxOutsideList(allIdx, mappedAugIdx)
		}else{
			chooseDip1 = mappedAugIdx[chooseDip1]
		}
		
		if(chooseDip2 == commonProbIdx){
			if(is.null(allIdx)) allIdx = 1:(2^snpLen)
			chooseDip2 = sampleIdxOutsideList(allIdx, mappedAugIdx)
		}else{
			chooseDip2 = mappedAugIdx[chooseDip2]
		}
		
	}
	
	reDip = range(chooseDip1, chooseDip2)
	
	return(reDip)
	
}


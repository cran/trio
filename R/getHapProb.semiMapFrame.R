getHapProb.semiMapFrame <-
function( semiMapFrame, resiProbCol, augIdxCol, probCol, snpLen){
	## SemiAugMap only have the major haplotypes' expression and prob,
	## only a number for prob for each of the other haplotypes 
	## SemiAugMap also has a column showing the index of the haplotype in the fixed AugMap
	ifD = F
	
	leftOverProb = semiMapFrame[1, resiProbCol]
	mapProb = semiMapFrame[ , probCol]
	mappedAugIdx = semiMapFrame[ , augIdxCol]
	
	commonProbIdx = length(mapProb)
	
	lef = leftOverProb/(2^snpLen-commonProbIdx)
	allProb = rep(lef, length=2^snpLen)
	
	allProb[mappedAugIdx]=mapProb
	
	return(allProb)
}


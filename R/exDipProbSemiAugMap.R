exDipProbSemiAugMap <-
function(dipMap, semiMapFrame, resiProbCol, augIdxCol, probCol, snpLen){
	## SemiAugMap only have the major haplotypes' expression and prob,
	## only a number for prob for each of the other haplotypes 
	## SemiAugMap also has a column showing the index of the haplotype in the fixed AugMap
	ifD = F
	
	dipProb = as.vector(dipMap)
	
	leftOverProb = semiMapFrame[1, resiProbCol]
	
	mappedAugIdx = semiMapFrame[ , augIdxCol]
	
	mapProb = semiMapFrame[ , probCol]
	
	commonProbIdx = length(mapProb) + 1
	
	curIdx = match(dipProb, mappedAugIdx, nomatch=commonProbIdx )
	fittedFilter = curIdx==commonProbIdx
	
	## find out how many un/specified haplotype will fit the observed data
	fittedInMap= unique(curIdx[!fittedFilter])
	fittedOutMap = unique(dipProb[fittedFilter])
	
	if(ifD) print("fittedInMap")
	if(ifD) print(fittedInMap)
	if(ifD) print("fittedOutMap")
	if(ifD) print(fittedOutMap)
	
	## restandarize the prob
	if(length(fittedOutMap)>0){
		mapProbRe=rep(0, commonProbIdx)
		mapProbRe[commonProbIdx]=leftOverProb/length(fittedOutMap)
		mapProbRe[fittedInMap] = mapProb[fittedInMap]/sum(mapProb[fittedInMap])*(1-leftOverProb)
		
	}else{
		mapProbRe =rep(0, commonProbIdx-1)
		mapProbRe[fittedInMap] = mapProb[fittedInMap]/sum(mapProb[fittedInMap])
		
	}
	
	
	if(ifD) print(mapProbRe)
	
	dipProb = mapProbRe[curIdx]
	
	if(ifD) print(dipProb)
	
	## remember the dipMap is passed in as a vector filted by column
	idxSeq = 1:(length(curIdx)/2)
	reDipProb = dipProb[idxSeq] * dipProb[idxSeq+ (length(curIdx)/2) ]
	
	return(reDipProb)
	
}


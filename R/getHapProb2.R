getHapProb2 <-
function(selIdx, semiMapFrame, resiProbCol, augIdxCol, probCol, snpLen, restandard=F){
	## SemiAugMap only have the major haplotypes' expression and prob,
	## only a number for prob for each of the other haplotypes 
	## SemiAugMap also has a column showing the index of the haplotype in the fixed AugMap
	ifD = F
	
	leftOverProb = semiMapFrame[1, resiProbCol]
	mapProb = semiMapFrame[ , probCol]
	mappedAugIdx = semiMapFrame[ , augIdxCol]
	
	commonProbIdx = length(mapProb)
	
	findMatchMajor = match(selIdx, mappedAugIdx)
	selMajor = findMatchMajor[!is.na(findMatchMajor)]
	if(ifD) print(semiMapFrame)
	if(ifD) print(selMajor)
	
	selMajor.ct = sum(!is.na(findMatchMajor))
	selMinor.ct = sum(is.na(findMatchMajor))
	
	if(commonProbIdx == (2^snpLen)){
		## if the map is complete
		
		if(selMinor.ct>0){
			stop (paste("non existing selIdx=[", paste(selIdx, collapse=";", sep=""), "]", sep=""))
		}else if(selMajor.ct>0){
			newProb = mapProb[ selMajor ]
		}else {
			stop (paste("not valide selIdx=[", paste(selIdx, collapse=";", sep=""), "]", sep=""))
		}
	}else{
		## if the map is incomplete
		
		augProb = c(mapProb, leftOverProb/(2^snpLen - commonProbIdx))
		idxMatched = selMajor
		if(selMinor.ct>0){
			validIdx = selIdx[is.na(findMatchMajor)]
			tt = (validIdx>=1 & validIdx <=2^snpLen)
			if(sum(tt)!=length(validIdx)){
				stop (paste("non existing selIdx=[", paste(selIdx, collapse=";", sep=""), "]", sep=""))
			}
			
			idxMatched = findMatchMajor
			idxMatched[is.na(findMatchMajor)] = commonProbIdx+1
			
		}
		newProb = augProb[idxMatched]
	}
	
	if(restandard){
		newProb = newProb/sum(newProb)
	}
	
	return(newProb)
}


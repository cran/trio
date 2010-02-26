sampleHapSemiAugMap2 <-
function(semiMapFrame, resiProbCol, augIdxCol, probCol, snpLen, exHapIdx=NULL){
	## SemiAugMap only have the major haplotypes' expression and prob,
	## only a number for prob for each of the other haplotypes 
	## SemiAugMap also has a column showing the index of the haplotype in the fixed AugMap
	ifD = F
	fStr = "[sampleHapSemiAugMap2]:"
	if(ifD) print(fStr)
	
	leftOverProb = semiMapFrame[1, resiProbCol]
	mapProb = semiMapFrame[ , probCol]
	mappedAugIdx = semiMapFrame[ , augIdxCol]
	
	exMajor = NULL
	exMinor = NULL
	exMajor.ct = 0
	exMinor.ct = 0
	commonProbIdx = length(mapProb)
	chooseHapIdx = NULL
	
	if(!is.null(exHapIdx)){
		
		matched = (  (exHapIdx >=1) & (exHapIdx <= (2^snpLen)))
		
		if(sum(matched)!=length(exHapIdx)) stop( paste("not valid exHapIdx=[", paste(exHapIdx, collapse=";", sep=""), "]", sep=""))
		
		majorMatch = is.element(exHapIdx, mappedAugIdx)
		exMajor = exHapIdx[majorMatch]
		exMinor = exHapIdx[!majorMatch]
		
		leftMajorRow = (1:commonProbIdx)[!is.element(mappedAugIdx, exMajor)]
		exMajor.ct = sum(majorMatch)
		exMinor.ct = sum(!majorMatch)
		
	}
	
	if(ifD){
		print(semiMapFrame)
		if(!is.null(exHapIdx)){
			print(paste("exMajor=[", paste(exMajor, collapse=";", sep=""), "]", sep=""))
			print(paste("leftMajorRow=[", paste(leftMajorRow, collapse=";", sep=""), "]", sep=""))
			print(paste("exMinor=[", paste(exMinor, collapse=";", sep=""), "]", sep=""))
			print(paste("majorMatch=[", paste(majorMatch, collapse=";", sep=""), "]", sep=""))
		}
	}
	
	
	## need to take account the situation where all the hap is presented.
	if(commonProbIdx == 2^snpLen){
		if(exMinor.ct>=1) stop(paste(fStr, " impossible to exclude minor when the map is completed."))
		## the only excluded are major
		if(exMajor.ct>0){
			## readjust the prob
			mapProbN = mapProb[leftMajorRow]
			mappedAugIdxN = mappedAugIdx[leftMajorRow]
			
			if(ifD) print(paste("mapProb=[", paste(round(mapProbN, 2), collapse=";", sep=""), "]", sep=""))
			if(ifD) print(paste("length=[", commonProbIdx-exMajor.ct , "]", sep=""))
			
			chooseRow = sample(1:(commonProbIdx-exMajor.ct), size=1, prob=mapProbN)
			chooseHapIdx = mappedAugIdxN[chooseRow]
			
		}else{
			## sample from completed map
			if(ifD) print(paste("mapProb=[", paste(round(mapProb, 2), collapse=";", sep=""), "]", sep=""))
			if(ifD) print(paste("length=[", commonProbIdx , "]", sep=""))
			
			chooseRow = sample(1:commonProbIdx, size=1, prob=mapProb)
			chooseHapIdx = mappedAugIdx[chooseRow]
		}
		
	}else{
		## not a complete map
		if(exMajor.ct==0 & exMinor.ct==0){
			## no exclusion
			
			mapProbN = c(mapProb,  leftOverProb)
			
			if(ifD) print("(exMajor.ct==0 & exMinor.ct==0):")
			if(ifD) print(paste("mapProbN=[", paste(round(mapProbN, 2), collapse=";", sep=""), "]", sep=""))
			if(ifD) print(paste("length=[", commonProbIdx+1 , "]", sep=""))
			
			chooseRow = sample(1:(commonProbIdx+1), size=1, prob=mapProbN)
			
			if(chooseRow==commonProbIdx+1){
				allIdx = 1:(2^snpLen)
				chooseHapIdx = sampleIdxOutsideList(allIdx, mappedAugIdx)
			}else{
				chooseHapIdx = mappedAugIdx[chooseRow]
			}
			
		}else if(exMajor.ct==0 & exMinor.ct>0){
			## only exclude minor
			before.minor = 2^snpLen - commonProbIdx
			mapProbN = c(mapProb,  leftOverProb/before.minor*(before.minor-exMinor.ct))
			
			if(ifD) print("(exMajor.ct==0 & exMinor.ct>0):")
			if(ifD) print(paste("mapProbN=[", paste(round(mapProbN, 2), collapse=";", sep=""), "]", sep=""))
			if(ifD) print(paste("length=[", commonProbIdx+1 , "]", sep=""))
			
			chooseRow = sample(1:(commonProbIdx+1), size=1, prob=mapProbN)
			
			
			if(chooseRow==commonProbIdx+1){
				allIdx = 1:(2^snpLen)
				chooseHapIdx = sampleIdxOutsideList(allIdx, c(mappedAugIdx, exMinor))
			}else{
				chooseHapIdx = mappedAugIdx[chooseRow]
			}      
			
		}else if(exMajor.ct>0 & exMinor.ct==0){
			## only exclude major
			mapProbN = mapProb[leftMajorRow]
			mappedAugIdxN = mappedAugIdx[leftMajorRow]
			mapProbN = c(mapProbN,  leftOverProb)
			
			if(ifD) print("(exMajor.ct>0 & exMinor.ct==0):")
			if(ifD) print(paste("mapProbN=[", paste(round(mapProbN, 2), collapse=";", sep=""), "]", sep=""))
			if(ifD) print(paste("length=[", commonProbIdx-exMajor.ct+1 , "]", sep=""))
			
			chooseRow = sample(1:(commonProbIdx-exMajor.ct+1), size=1, prob=mapProbN)
			
			if(chooseRow==(commonProbIdx-exMajor.ct+1)){
				allIdx = 1:(2^snpLen)
				## tricky, can only sample minor, augIdx pool are not reduced
				chooseHapIdx = sampleIdxOutsideList(allIdx, c(mappedAugIdx))
			}else{
				chooseHapIdx = mappedAugIdxN[chooseRow]
			}      
		}else if(exMajor.ct>0 & exMinor.ct>0){
			## exlude both major and minor
			mapProbN = mapProb[leftMajorRow]
			mappedAugIdxN = mappedAugIdx[leftMajorRow]
			
			before.minor = 2^snpLen - commonProbIdx
			
			mapProbN = c(mapProbN,  leftOverProb/before.minor*(before.minor-exMinor.ct))
			
			
			if(ifD) print("(exMajor.ct>0 & exMinor.ct>0):")
			if(ifD) print(paste("mapProbN=[", paste(round(mapProbN, 2), collapse=";", sep=""), "]", sep=""))
			if(ifD) print(paste("length=[", commonProbIdx-exMajor.ct+1 , "]", sep=""))
			
			chooseRow = sample(1:(commonProbIdx-exMajor.ct+1), size=1, prob=mapProbN)
			
			if(chooseRow==(commonProbIdx-exMajor.ct+1)){
				allIdx = 1:(2^snpLen)
				## tricky, can only sample minor, augIdx pool are not reduced
				chooseHapIdx = sampleIdxOutsideList(allIdx, c(mappedAugIdx, exMinor))
			}else{
				chooseHapIdx = mappedAugIdxN[chooseRow]
			}    
			
			
		}else{
			stop("Nothing here")
		}
		
	}
	
	return(chooseHapIdx)
	
}


fHapBkIdx2DipTb <-
function(appVarNames, filteredBkInfo, idxList, snpCoding, reqIn=NULL, reqDigits = NULL, expression, snpLen=nchar(expression)){
	ifD = F
	fStr = "[fHapBkIdx2DipTb]:"
	if(ifD) {
		print(qp(fStr, "start"))
		print(filteredBkInfo)
	}
	if( min( snpCoding==c(0,1,2,3)) <1 )  
		stop (paste("\nData configuration is not right:\n", "snpCoding=[", paste(snpCoding, collapse=";", sep=""), "]", sep=""))
	
	
	hetoIn = filteredBkInfo$hetoIn
	
	if(!is.null(hetoIn)){
		## need to build up the dip pair following the heto digit requirment\
		bkCt = length(idxList)
		
		if(bkCt<2)
			stop(paste("\nCannot find more than two haplotype to build heto digit for data:(", expression, ").", sep=""))
		hetoSeq = hetoIn
	}else{
		## only need to process the digit requirement posed by child
		hetoSeq = NULL
		bkCt = length(idxList)
	}
	
	if(ifD & (!is.null(hetoSeq))) print(paste("hetoSeq:", paste(hetoSeq, collapse=";", sep="")))
	
	## if no heto digit and no requirment posed by child, just return the all possible dip pair
	if(is.null(reqIn) & is.null(hetoIn) ){
		dipMap = util.it.upTriCombIdx(idxList, diag=T, re.ordered=T)
		#print(dipMap)
		return(dipMap)
	}
	
	
	## if no heto digit but some requirement posed by child
	if(is.null(hetoIn) & (!is.null(reqIn)) ){
		dipMap = NULL
		for ( i in 1:bkCt ){
			curBkIdx = idxList[i]
			## because no heto digits restriction, then dip built by the same hap should be considered.
			j = i
			
			while ( j <= bkCt ){
				## compare the heto digits
				othBkIdx = idxList[j]
				
				## check the required digits
				meetReq = T
				tmpReq = 1
				while( tmpReq <= length(reqIn)){
					curDigitReq = reqDigits[tmpReq]
					tmpSeq = seq.int(from=1, to=2^(reqIn[tmpReq]-1), by=1)
					oneMeetReq = isDigitAtLociIdx(hapIdx=curBkIdx, digit=curDigitReq,
							lociCt=snpLen, lociIdx=reqIn[tmpReq], intVec = tmpSeq)
					othMeetReq = isDigitAtLociIdx(hapIdx=othBkIdx, digit=curDigitReq,
							lociCt=snpLen, lociIdx=reqIn[tmpReq], intVec = tmpSeq)
					
					meetReq = meetReq & (oneMeetReq | othMeetReq)
					tmpReq = tmpReq + 1
					if(!meetReq) tmpReq = length(reqIn)+10
					
				}
				## keep the matched pair
				if(meetReq)  dipMap = rbind(dipMap, range(curBkIdx, othBkIdx))
				j = j + 1
			} # while ( j <= bkCt ){
		} # for ( i in 1:bkCt ){
		if(is.null(dipMap))
			stop(paste("\nCannot find the haplotype pairs for parent meet digit requirement, exp=(", expression, ").", sep=""))
		#print(dipMap)
		return(dipMap)   
	}
	
	## if( !is.null(hetoIn) & [ is.null(reqIn) | !is.null(reqIn) ] )
	
	idx4hapDigit = NULL
	## check out the global variables
	tryCatch({
				
				tmpGetObj = NULL
				tmpGetObj = get(appVarNames$digit, env=.GlobalEnv)
				
				idx4hapDigit$digitMap1 = tmpGetObj$digitMap1[1:(2^(snpLen-1)), 1:snpLen]
				idx4hapDigit$digitMap2 = tmpGetObj$digitMap2[1:(2^(snpLen-1)), 1:snpLen]
				
				rm(tmpGetObj)
				##gc()
			}, error=function(e){
				errTrace = paste(e, collapse=";", sep="") 
				stop(paste("\n", fStr, errTrace, "\nApp-wise Global Variable ", appVarNames$digit, " does not exisit."))
			})
	
	
	dipMap = NULL
	## check for all required heto
	
	## all the idxList meet the homo digit requirement, but for heto digits
	## we need to match up the pair for all the heto digits
	#print(idxList)
	#print(idx4hapDigit)
	for ( i in 1:bkCt ){
		curBkIdx = idxList[i]
		## because heto digits restriction, then dip built by the same hap should be excluded.
		j = i + 1 
		while( j <= bkCt ){
			## compare the heto digits
			othBkIdx = idxList[j]
			
			if(ifD) print(paste("i=", i, "; j=", j, sep=""))
			
			## suppose cur bk contribute a digit 1
			curLociMatch1 = util.matrix.colIdx4Match(ma=idx4hapDigit$digitMap1[, hetoSeq, drop=F], val=curBkIdx)
			## suppose oth bk contribute a digit 2
			othLociMatch2 = util.matrix.colIdx4Match(ma=idx4hapDigit$digitMap2[, hetoSeq, drop=F], val=othBkIdx)
			
			## suppose oth bk contribute a digit 1
			othLociMatch1 = util.matrix.colIdx4Match(ma=idx4hapDigit$digitMap1[, hetoSeq, drop=F], val=othBkIdx)
			## suppose cur bk contribute a digit 2
			curLociMatch2 = util.matrix.colIdx4Match(ma=idx4hapDigit$digitMap2[, hetoSeq, drop=F], val=curBkIdx)               
			
			#print(curLociMatch1)
			#print(othLociMatch2)
			
			#print(othLociMatch1)
			#print(curLociMatch2)
			
			## change code, may cause trouble      
			##        if(length(curLociMatch1)==0) curLociMatch1 = 100
			##        if(length(curLociMatch2)==0) curLociMatch2 = 101
			##        if(length(othLociMatch1)==0) othLociMatch1 = 102
			##        if(length(othLociMatch2)==0) othLociMatch2 = 103
			##        tmp = sum(c( suppressWarnings(curLociMatch1==othLociMatch2),
			##                    suppressWarnings(curLociMatch2==othLociMatch1)))
			
			
			## check same length
			tmp=0
			
			if(length(curLociMatch1)==length(othLociMatch2) ){
				if(length(curLociMatch1)!=0 ){
					tmp=sum(curLociMatch1==othLociMatch2)
				}
			}
			
			if(length(curLociMatch2)==length(othLociMatch1) ){
				if(length(curLociMatch2)!=0 ){
					tmp=tmp+sum(curLociMatch2==othLociMatch1)
				}
			}
			
			
			if(tmp!=0){
				
				if(tmp == length(hetoSeq)){
					## find the matched pair
					## check the required digit
					if(is.null(reqIn)){
						## keep the matched pair
						dipMap = rbind(dipMap, range(curBkIdx, othBkIdx))
						if(ifD) print(paste("curBkIdx=", curBkIdx, ": othBkIdx=", othBkIdx, sep=""))
					}else{
						## check the required digits
						meetReq = T
						tmpReq = 1
						while( tmpReq <= length(reqIn)){
							curDigitReq = reqDigits[tmpReq]
							tmpSeq = seq.int(from=1, to=2^(reqIn[tmpReq]-1), by=1)
							oneMeetReq = isDigitAtLociIdx(hapIdx=curBkIdx, digit=curDigitReq,
									lociCt=snpLen, lociIdx=reqIn[tmpReq], intVec = tmpSeq)
							othMeetReq = isDigitAtLociIdx(hapIdx=othBkIdx, digit=curDigitReq,
									lociCt=snpLen, lociIdx=reqIn[tmpReq], intVec = tmpSeq)
							
							meetReq = meetReq & (oneMeetReq | othMeetReq)
							tmpReq = tmpReq + 1
							if(meetReq) tmpReq = length(reqIn)+10
							
						}
						## keep the matched pair
						if(meetReq)  {
							dipMap = rbind(dipMap, range(curBkIdx, othBkIdx))
							if(ifD) print(paste("curBkIdx=", curBkIdx, "; othBkIdx=", othBkIdx, sep=""))
						}
					} # if(is.null(reqIn)){
				} # if(tmp == length(hetoSeq)){
			} #  if(tmp!=0){
			
			j = j + 1
			
		} # while( j <= bkCt ){
	} # for ( i in 1:bkCt ){
	rm(idx4hapDigit)
	##gc()
	if(is.null(dipMap))
		stop(paste("\n Cannot find the haplotype pairs for parent meet heto/digit requirement, exp=(", expression, ").", sep=""))
	
	#print(dipMap)
	
	return(dipMap)
	
}


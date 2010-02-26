procSemiAugMap <-
function(appVarNames, homoHetoInfo, snpLen){
	
	fStr="[procSemiAugMap]:"
	
	if(is.null(homoHetoInfo$homoIn)){
		## if not homo digit requirment, every hap is possible
		retainList = 1:(2^snpLen)
		return(retainList)
	}
	
	## if there is homo digit requirement
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
				## traceback()
				## print(e)
				errTrace = paste(e, collapse=";", sep="")
				stop(paste("\n", fStr, errTrace, "\nApp-wise Global Variable ", appVarNames$digit, " does not exisit."))
			})
	
	
	
	## first process the homo digit
	
	retainList = NULL
	matchedIdx = NULL
	for( i in 1:length(homoHetoInfo$homoIn)){
		curDigit =  homoHetoInfo$homoDigit[i]
		tmpPos = homoHetoInfo$homoIn[i]
		if(curDigit == 1) {
			matchedIdx = idx4hapDigit$digitMap1[, tmpPos, drop=F]
			if(is.null(retainList)){
				retainList = matchedIdx
			}else{
				retainList = retainList[is.element(retainList, matchedIdx)]
			}
			## print(retainList)
		}else{
			matchedIdx = idx4hapDigit$digitMap2[, tmpPos, drop=F]
			if(is.null(retainList)){
				retainList = matchedIdx
			}else{
				retainList = retainList[is.element(retainList, matchedIdx)]
			}
			##  print(retainList)
		} ## if(curDigit == 1) {
	} ## for( i in homoHetoInfo$homoIn){
	
	rm(idx4hapDigit)
	##gc()
	if( is.null (retainList) )
		stop(paste("\nBlock infomation error. No population haplotype matches existing data:(", expression, ").", sep=""))
	
	return (retainList)
	
}


ESp.impu1Par.C.sp <-
function(hap, hap.p, straSeq, semiMapFrame, resiProbCol, augIdxCol, probCol, snpLen){
	#print(paste("straSeq=", straSeq))
	## give out 6 item vectors for the hap pairs for imputed par, given par and child
	if(straSeq==1){
		## i,i|ij|ii
		re = c(hap, hap,   sort(c(hap, hap.p)), hap, hap)
	}
	
	if(straSeq==2){
		## i, !=ij|ij|ii
		sp = sampleHapSemiAugMap2(semiMapFrame, resiProbCol, augIdxCol, probCol, snpLen, exHapIdx=c(hap, hap.p))
		#print(paste("sp=", sp))
		re = c(sort(c(hap, sp)), sort(c(hap, hap.p)), hap, hap)
	}
	if(straSeq==3){
		## ij|ij|ii
		re = c(sort(c(hap, hap.p)), sort(c(hap, hap.p)), hap, hap)
	}
	
	return(re)
	
}


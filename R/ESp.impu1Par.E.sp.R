ESp.impu1Par.E.sp <-
function(hap, hap.p, straSeq, semiMapFrame, resiProbCol, augIdxCol, probCol, snpLen ){
	# print(straSeq)
	re = NULL
	## give out 6 item vectors for the hap pairs for imputed par, given par and child
	if(straSeq==1){
		## k,k|ij|ik
		re = c(hap[2], hap[2], sort(c(hap[1], hap.p)), sort(hap))
	}
	
	if(straSeq==2){
		## k, !=k|ij|ik
		sp = sampleHapSemiAugMap2(semiMapFrame, resiProbCol, augIdxCol, probCol, snpLen, exHapIdx=c(hap[2]))
		re = c(sort(c(hap[2], sp)),     sort(c(hap[1], hap.p)), sort(hap))
	}
	
	
	return(re)
	
}


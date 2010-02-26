ESp.impu1Par.D.sp <-
function(hap,  straSeq, semiMapFrame, resiProbCol, augIdxCol, probCol, snpLen ){
	
	## give out 6 item vectors for the hap pairs for imputed par, given par and child
	if(straSeq==1){
		## kk|ii|ik
		re = c(hap[2], hap[2], hap[1], hap[1], sort(hap))
	}
	
	if(straSeq==2){
		## k, !=k|ii|ik
		sp = sampleHapSemiAugMap2(semiMapFrame, resiProbCol, augIdxCol, probCol, snpLen, exHapIdx=c(hap[2]))
		re = c(sort(c(hap[2], sp)),     hap[1], hap[1], sort(hap))
	}
	
	return(re)
	
}


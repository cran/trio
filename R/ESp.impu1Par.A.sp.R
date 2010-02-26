ESp.impu1Par.A.sp <-
function(hap,  straSeq, semiMapFrame, resiProbCol, augIdxCol, probCol, snpLen ){
	
	## give out 6 item vectors for the hap pairs for imputed par, given par and child
	if(straSeq==1){
		## i,i|ii|ii
		re = rep(hap, 6)
	}
	
	if(straSeq==2){
		## i, !=i|ii|ii
		sp = sampleHapSemiAugMap2(semiMapFrame, resiProbCol, augIdxCol, probCol, snpLen, exHapIdx=hap)
		re = c(sort(c(hap, sp)), rep(sort(hap), 4))
	}
	
	return(re)
	
}


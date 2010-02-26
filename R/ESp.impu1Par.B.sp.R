ESp.impu1Par.B.sp <-
function(hap,  straSeq, semiMapFrame, resiProbCol, augIdxCol, probCol, snpLen ){
	
	## give out 6 item vectors for the hap pairs for imputed par, given par and child
	if(straSeq==1){
		## i,i|ik|ik
		re = c(rep(hap[1], 2), sort(hap), sort(hap))
	}
	
	if(straSeq==2){
		## i, !=ik|ik|ik
		sp = sampleHapSemiAugMap2(semiMapFrame, resiProbCol, augIdxCol, probCol, snpLen, exHapIdx=hap)
		re = c(sort(c(hap[1], sp)), sort(hap), sort(hap))
	}
	if(straSeq==3){
		## ik|ik|ik
		re = rep(sort(hap), 3)
	}
	
	if(straSeq==4){
		## k, !=ik|ik|ik
		sp =  sampleHapSemiAugMap2(semiMapFrame, resiProbCol, augIdxCol, probCol, snpLen, exHapIdx=hap)    
		re = c(sort(c(hap[2], sp)),  sort(hap), sort(hap))   
	}
	if(straSeq==5){
		## k,k|ik|ik
		re = c(rep(hap[2], 2),  sort(hap), sort(hap))
	}
	
	return(re)
	
}


HRCB.Esp1Rule.sampleKid <-
function(rule, caseNo, spStrata, supHapProb){
	ifD =F
	fn = "HRCB.Esp1Rule.sampleKid::"
	
	tmpObj = NULL
	if( is.character(spStrata)){
		
		a = load(paste(spStrata, ".RData", sep=""))
		tmpObj = get(a[1])
		
	}else{
		tmpObj = spStrata
	}
	
	
	HRCBStra.A = tmpObj$HRCBStra.A
	HRCBStra.B = tmpObj$HRCBStra.B
	HRCBStra.C = tmpObj$HRCBStra.C
	HRCBStra.D = tmpObj$HRCBStra.D
	
	signalRule = rule$slist[[1]]
	inter = rule$inter
	coef = signalRule$coef
	
	HR = exp(inter+coef)/(1+exp(inter+coef))
	LR = exp(inter)/(1+exp(inter))
	
	
	if(HRCBStra.A$ct>0) {
		risk.A = HRCBStra.A$idProb[,2]*unlist(HR)
	}else{
		risk.A = 0
	}
	if(HRCBStra.B$ct>0) {
		risk.B = HRCBStra.B$idProb[,2]*unlist(HR)
	}else{
		risk.B = 0
	}
	if(HRCBStra.C$ct>0) {
		risk.C = HRCBStra.C$idProb[,2]*unlist(LR)
	}else{
		risk.C = 0
	}
	risk.D = HRCBStra.D$idProb[,2]*unlist(LR)
	
	risk.sumHR = lapply(list(risk.A, risk.B), FUN=sum)
	risk.sumLR = lapply(list(risk.C, risk.D), FUN=sum)
	
	
	# check the sum of prob
	#print(paste(fn, " population risk =", sum(unlist(c(risk.sumHR, risk.sumLR)))))
	
	
	stra.sampled = rmultinom(1, size=caseNo, prob=c(risk.sumHR, risk.sumLR) )
	stra.cumsam = cumsum(stra.sampled)  
	if(ifD){
		print(paste("sampled group:", paste(stra.sampled, collapse="; ")))
	}
	
	
	simMa = matrix(NA, nrow=caseNo, ncol=2)
	## simu for group A
	if(stra.sampled[1]>0){
		#print("a")
		simMa[1:stra.cumsam[1],] = HRCBSpGrp.sp(grp=HRCBStra.A, size=stra.sampled[1], hapProb=supHapProb, riskProb=risk.A)
		#print("a")
	}
	## simu for group B
	if(stra.sampled[2]>0){
		#print("b")
		simMa[(stra.cumsam[1]+1):stra.cumsam[2] ,] = HRCBSpGrp.sp(grp=HRCBStra.B, size=stra.sampled[2], hapProb=supHapProb, riskProb=risk.B)
		
	}
	## simu for group C
	if(stra.sampled[3]>0){
		#print("c")
		simMa[(stra.cumsam[2]+1):stra.cumsam[3], ] =
				HRCBSpGrp.sp(grp=HRCBStra.C, size=stra.sampled[3], hapProb=supHapProb)
		
	}
	## simu for group D
	if(stra.sampled[4]>0){
		#print("d")
		simMa[(stra.cumsam[3]+1):stra.cumsam[4], ] =
				HRCBSpGrp.sp(grp=HRCBStra.D, size=stra.sampled[4], hapProb=supHapProb, grpA =HRCBStra.A )
		
	}
	
	return(simMa)
}


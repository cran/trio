ESp.impu1Par.E <-
function(hap, prob.p, semiMapFrame, resiProbCol, augIdxCol, probCol, snpLen ){
	
	prob.c = getHapProb2(selIdx=hap, semiMapFrame,
			resiProbCol=resiProbCol,
			augIdxCol=augIdxCol,
			probCol=probCol,
			snpLen, restandard=F)
	
	
	stra.1 = (prob.c^2)*(prob.p[1]*prob.p[2])
	stra.2 = (prob.c)*(1-prob.c)*(prob.p[1]*prob.p[2])
	
	return(c(stra.1, stra.2))
	
}


ESp.impu1Par.B <-
function(hap, prob.p, semiMapFrame, resiProbCol, augIdxCol, probCol, snpLen){
	
	prob.c1 = getHapProb2(selIdx=hap[1], semiMapFrame,
			resiProbCol=resiProbCol,
			augIdxCol=augIdxCol,
			probCol=probCol,
			snpLen, restandard=F)
	prob.c2 = getHapProb2(selIdx=hap[2], semiMapFrame,
			resiProbCol=resiProbCol,
			augIdxCol=augIdxCol,
			probCol=probCol,
			snpLen, restandard=F)
	
	stra.1 = (prob.c1^2)*(prob.p[1]*prob.p[2])
	stra.2 = (prob.c1)*(1-prob.c1-prob.c2)*(prob.p[1]*prob.p[2])
	stra.3 = 2*(prob.c1*prob.c2)*(prob.p[1]*prob.p[2])
	stra.4 = (prob.c2)*(1-prob.c1-prob.c2)*(prob.p[1]*prob.p[2])
	stra.5 = (prob.c2^2)*(prob.p[1]*prob.p[2])
	
	return(c(stra.1, stra.2, stra.3, stra.4, stra.5))
	
}


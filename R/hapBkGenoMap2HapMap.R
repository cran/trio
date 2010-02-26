hapBkGenoMap2HapMap <-
function(hapBkGenoMap, reSimuMap=T){
	
	singleton = hapBkGenoMap$genomeMarkerInfo[hapBkGenoMap$genomeMarkerInfo$qu.type==1, 2 ]
	singletonGeno = hapBkGenoMap$genoOnlyMap$bks[ singleton ]
	newFrame = NULL
	for( i in 1: (length(singletonGeno))){
		g = singletonGeno[[i]]
		newF = geno2hapFreq(prekey=g[1,hapBkGenoMap$genoOnlyMap$chCol] ,
				block =g[1,hapBkGenoMap$genoOnlyMap$blockCol] ,
				genoFreq=g[,hapBkGenoMap$genoOnlyMap$probCol] ,
				snpBase=hapBkGenoMap$hapBkOnlyMap$snpBase,
				snpSeq=singleton[i])
		newFrame=c(newFrame, list(newF))
	}
	
	singletonSeqId = which(hapBkGenoMap$genomeMarkerInfo$qu.type==1)
	allSeqId = 1: (nrow(hapBkGenoMap$genomeMarkerInfo))
	leftSeqId = allSeqId[ -singletonSeqId ]
	
	## bind with old
	allBkMap = c(hapBkGenoMap$hapBkOnlyMap$bks, newFrame)
	## reshuffle the bks
	newBkMap = allBkMap[ match(allSeqId, c(leftSeqId, singletonSeqId)) ]
	
	
	## get new DF
	newDf = matrix(NA, nrow=length(singletonGeno)*2+nrow(hapBkGenoMap$hapBkOnlyMap$df), ncol=7)
	newDf = data.frame(key=I(rep("-", nrow(newDf))), newDf)
	
	rowIdx = 0
	for( j in 1:(length(newBkMap))){
		hapDf = newBkMap[[j]]
		newDf[ (rowIdx+1):(rowIdx+nrow(hapDf)), ]= hapDf
		rowIdx = rowIdx+nrow(hapDf)
	}
	
	#str(newDf)
	if(reSimuMap){
		simuDf = newDf[ , c(1, 4, 5)]
		simuMap = bkMap.constr(data=simuDf, keyCol=1, hapLenCol=NULL, expCol=2, probCol=3)
		return(list( newBkMap=newBkMap, newDf=newDf, simuMap=simuMap))
	}
	
	return(list(  newBkMap=newBkMap, newDf=newDf ))
}


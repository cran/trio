bkMap.ESp.apply1Rule.stepBy <-
function(bkMap, signalRule){
	ifD = F
	fN = "bkMap.ESp.apply1Rule"
	if(ifD) print(paste(fN, "begin::"))
	
	signal = signalRule$signal
	
	total.bkct = length( bkMap$snpCt )
	
	# construct the beginning/ending index of snp in original map
	idx.be = matrix(NA, nrow=total.bkct, ncol=2)
	idx.be[,2]=cumsum(bkMap$snpCt)
	idx.be[,1]=c(1, cumsum(bkMap$snpCt)+1)[1:total.bkct]
	
	
	causalBkUniqOrderedIdx = bkMap.findHRCBIdx(bkMap, as.integer(signalRule$var.idx), re.keys=F, unique=T)
	## find out # of haplotype in each HRCB
	HRCB.hapCt = bkMap$bkLens [causalBkUniqOrderedIdx]
	HRCB.hapCtProd = cumprod(HRCB.hapCt)[length(causalBkUniqOrderedIdx)]
	
	causalBkIdx.ruleOrder = bkMap.findHRCBIdx(bkMap, as.integer(signalRule$var.idx), re.keys=F, unique=F)
	
	#if(length(causalBkUniqOrderedIdx)==1){
	
	if(length(causalBkIdx.ruleOrder)==1){
		
	}else{
		## need to update the rule a little
		if(ifD) print( binaTree.toStr(signal))
		signal.n = binaTree.shuffleNode(signal)    
		if(ifD) print( binaTree.toStr(signal.n))
		signal.f = binaTree.merge(signal.n)
		if(ifD) print( binaTree.toStr(signal.f))
		
	}
	
	## have to find the matching with the HRCB only bkMap
	#causalBkIdx.ruleOrder = match(causalBkIdx.ruleOrder, causalBkIdx)
	
	elmMList = list()
	# for each block, need to find the hap pairs that are associated with higher risk
	
	for ( j in 1:length(causalBkIdx.ruleOrder)){
		bkIdx = causalBkIdx.ruleOrder[ j ]
		if(ifD) print(paste("bkIdx inside the bkMap =", bkIdx ))
		
		##  20Dec08Change!!!: straighten coding: 1, 2, 3 for 1-digitCoding, and 1, 2 for 2-digit
		
		genoMap =  hap2GenoBlock(hapExp = as.character(bkMap$bks[[bkIdx]][, bkMap$expCol]),
				hapProb = bkMap$bks[[bkIdx]][, bkMap$probCol], snpCt = bkMap$snpCt[bkIdx],
				alleleCode = bkMap$alleleCode )
		if(ifD) print(genoMap)
		item = idx.be[bkIdx,]
		## use recycling rules
		re = paste("v", rep(item[1]:item[2], each=2), sep="")
		newData.col = paste(re,  c("a", "b"), sep="")
		
		## singleSNPRule follow the original order in the genome
		## HARD CODE!!!HARD CODE: genoMa from hap2GenoBlock keep digits  1, 2 for 2-digit, need to change to 0/1 binary coding
		## 20Dec08Change!!! to use 2- for 0/1 binay coding, so the first alleleCode (changed to 1) correspond to the minor allele
		## and the second  (changed to 2) for major allele. Just need to coordinate with signal Rule.
		
		varIdx = match (signalRule$signal$elm.vlist[[j]], newData.col)
		treePred = (2-genoMap$genoMa)[,varIdx]
		## if the variable is not proceed by "not" operator
		## CHANGED on 8JAN09, from ==0 to !=-1. -1 means not
		if (signalRule$signal$elm.vMa[j,3]!=-1){
			nct = sum(treePred==1)
			if(ifD) print("not proceed by not")
			pairHRCB = genoMap$tb[  treePred==1, c(1,2), drop=F ]
		}else{
			## if the variable is proceed by "not" operator
			nct = sum(treePred==0)
			if(ifD) print("proceed by not")
			pairHRCB = genoMap$tb[  treePred==0, c(1,2), drop=F ]
		}
		
		#print(paste("matching hap pairs in block ", bkIdx))
		#print(pairHRCB)
		
		
		elmMList=c(elmMList, list(pairHRCB))
	} # for ( j in 1:length(causalBkIdx.new)){
	#print("This is elmMList")
	#print(elmMList)
	#print(causalBkIdx.ruleOrder)
	#print("qing mark")
	
	## simplied if only one HRCB
	
	if (length(causalBkIdx.ruleOrder)==1){
		
		superDipIdx = apply(pairHRCB, 1, FUN=util.it.triMatch, len=HRCB.hapCtProd)
		
		superDipIdx = sort(superDipIdx)
		
		## FIX LATER!!!FIX : can further reduce the workload
		HRCBGrp = t(sapply(superDipIdx, FUN=util.it.triMatch2, len=HRCB.hapCtProd, re.homo=T))
		
		# print(str(HRCBGrp))
		
		HRCBGrp.A = HRCBGrp[HRCBGrp[,3]==1, c(1,2), drop=F]
		
		HRCBGrp.BIdx = HRCBGrp[HRCBGrp[,3]==0,1]
		
		#print(HRCBGrp.A)
		#print(HRCBGrp.BIdx)
		return(list(A=HRCBGrp.A, B=HRCBGrp.BIdx))
		
	}
	
	hapGrids = binaTree.patternForm(binaTree=signal.f, elmMList=elmMList, bkIdx.RuleOrder=causalBkIdx.ruleOrder)
	if(ifD) print(str(hapGrids))
	##write.csv(hapGrids, file="hapGrids.csv") 
	
	## check whether no matching BK
	apply(hapGrids, 2, FUN=function(coo){
				qing.check = mean(is.na(coo))
				if(qing.check==1) stop(qp(fN, ": no matching pattern for one block"))
			})
	
	
	## need to find super-hap index for the matching patterns
	uniBkCt = ncol(hapGrids) / 2
	
	HRCB.hapCt2 = HRCB.hapCt
	
	if(ifD) print(hapGrids)
	supIdx = NULL
	
	for(n in 1:nrow(hapGrids)){
		
		if(ifD) {
			print(paste("proc row for pattern: row=", nrow(hapGrids)))
			#print(n)
			#print(hapGrids[n,])
		}
		if( n>240) {
			#print(hapGrids[n,])
		}
		hap1 = hapGrids[n, 1:uniBkCt] 
		m1 = filterHaps2SHaps(hap1, HRCB.hapCt2) 
		
		hap2 = hapGrids[n, (uniBkCt+1):(2*uniBkCt)]
		m2 = filterHaps2SHaps(hap2, HRCB.hapCt2)  
		
		hapPairs = qExpandTable(listOfFactor =list(m1, m2), removedRowIdx=NULL, re.row=F)
		hapPairs = t(apply(hapPairs, 1, FUN=range))
		superDipIdx = apply(hapPairs, 1, FUN=util.it.triMatch, len=HRCB.hapCtProd)
		superDipIdx = unique(unlist(superDipIdx))
		
		supIdx = c(supIdx, superDipIdx)
		supIdx = unique(supIdx)
		
		if(length(supIdx)>10^10) stop( qp(fN,": Too many matching super-haplotypes, exceeding 10^10."))
		#if(ifD) print(supIdx)
	}
	
	supDipAll = NULL
	if( HRCB.hapCtProd < 250 ){
		supDipAll = supIdx[sort.list(as.integer(supIdx), method="radix")]
	}else{
		supDipAll = supIdx[sort.list(as.integer(supIdx), method="quick", na.last=NA)]
	}
	
	
	#print(paste(fN, "QingMark2::sorted superDip:"))
	#print(length(supDipAll))
	## parse out the four stratum
	if (length(supDipAll)==0){
		warning(paste(fN, ":no matching allele on the HRCB"))
	}
	
	#print(supDipAll)
	## FIX LATER!!!FIX : can further reduce the workload
	HRCBGrp = t(sapply(supDipAll, FUN=util.it.triMatch2, len=HRCB.hapCtProd, re.homo=T))
	
	
	#print(HRCBGrp)
	# print(str(HRCBGrp))
	
	HRCBGrp.A = HRCBGrp[HRCBGrp[,3]==1, c(1,2), drop=F]
	
	HRCBGrp.BIdx = HRCBGrp[HRCBGrp[,3]==0,1]
	#save(HRCBGrp.A, file="HRCBGrp.A.RData")
	#save(HRCBGrp.BIdx, file="HRCBGrp.BIdx.RData")
	
	return(list(A=HRCBGrp.A, B=HRCBGrp.BIdx))
	
}


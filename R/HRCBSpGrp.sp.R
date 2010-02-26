HRCBSpGrp.sp <-
function(grp, size, hapProb, riskProb=NULL, grpA=NULL){
	ifD = F
	samMa = matrix(NA, nrow=size, ncol=2)
	#print(grp)
	#print(hapProb)
	
	
	straCt = nrow(grp$idProb)
	stra.sampled=NULL
	if(grp$type=="A" | grp$type=="B"){
		if (is.null(riskProb))   riskProb=grp$idProb[,2]
		stra.sampled = rmultinom(1, size=size, prob= riskProb )
	}else{
		stra.sampled = rmultinom(1, size=size, prob= grp$idProb[,2])
	}
	if(ifD) print(grp$type)
	if(ifD) print("strata ct")
	if(ifD) print(stra.sampled)
	
	#print("sampled")
	tmp.ma = cbind(1:straCt, stra.sampled)
	tmp.ma = tmp.ma[tmp.ma[,2]>0, ,drop=F]
	#print(tmp.ma)
	
	samMa[,1] =  unlist(apply(tmp.ma, 1, FUN = function(row, ma){
						rep(ma[row[1]], times=row[2])}, ma=grp$idProb[,1]))
	
	#print(samMa)
	
	if (grp$type=="A"){
		#id2 = grp$id2
		samMa[,2] = unlist(apply(tmp.ma, 1, FUN=function(row, id2, hapProb, idxct ){
							#print(hapProb)
							if(ifD) print(row)
							if(row[1]!=length(idxct)){
								tmp.range = idxct[row[1]]:(idxct[row[1]+1]-1)
							}else{
								tmp.range = idxct[row[1]]:length(id2)
							}
							tmp.sam = id2 [tmp.range]
							
							id2.sel = resample(tmp.sam, size=row[2], prob = hapProb[tmp.sam], replace=T)
							#print(id2.sel)
							if(length(id2.sel)==1) id2.sel = rep(id2.sel, times=row[2])
							#print(id2.sel)
							return(id2.sel)
							
						},  id2 = grp$id2, hapProb=hapProb, idxct = grp$idProb[,3]))
		if(ifD) print(samMa)
		samMa = cbind(pmin(samMa[,1], samMa[,2]), pmax(samMa[,1], samMa[,2]))
	}
	if (grp$type=="B" | grp$type=="C"){
		samMa[,2] = samMa[,1]
	}
	
	if (grp$type=="D"){
		## sample the second idx by rejection
		if(is.null(grpA)) stop("HRCBSpGrp.sam::error")
		samMa[,2] = unlist(
				apply(tmp.ma, 1, FUN = function(row, grpAid, idxct, id2){
							maByIdx = match(row[1], grpAid)
							if( !is.na(maByIdx)){
								if(ifD) print(row)
								if(maByIdx!=length(idxct)){
									tmp.range = idxct[maByIdx]:(idxct[maByIdx+1]-1)
								}else{
									tmp.range = idxct[maByIdx]:length(id2)
								}
								## ask the program to trace up
								tmp.sam = id2 [tmp.range]
								removeList = c(tmp.sam, row[1])
								
								# sample for the entire list, except the one already in A
								idx2 = sam.reject(allIdx=grp$idProb[,1], idxProb = hapProb, rejectIdx = removeList, size=row[2])
							}else{
								# sample for the others, exept the one that is the same as the column 1           
								idx2 = sam.reject(allIdx=grp$idProb[,1], idxProb = hapProb, rejectIdx = row[1], size=row[2])
							}
							return(idx2)
						}, grpAid = grpA$idProb[,1], idxct=grpA$idProb[,3], id2=grpA$id2))
		
		if(ifD) print(samMa)
		samMa = cbind(pmin(samMa[,1], samMa[,2]), pmax(samMa[,1], samMa[,2]))
		if(ifD) print(samMa)
	}  
	
	return(samMa)
	
}


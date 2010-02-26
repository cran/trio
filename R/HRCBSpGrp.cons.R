HRCBSpGrp.cons <-
function(hapProb, ids, type="A"){
	ifD = F
	## contain the first index and the haplotype probabilities for hap1, hap2 and joint
	idProb = NULL
	id2 = NULL
	
	if (type == "A"){
		## for type A, the ids are two column indexes
		## need to find the list of second indexes for the same first index
		tmp.row = nrow(ids)
		if (tmp.row==0) {
			grpA = list(type="A",  ct=0)
			
		}else{
			
			idsAll = matrix(NA, nrow=tmp.row*2, ncol=2)
			idsAll[1:tmp.row,]=ids
			idsAll[(1:tmp.row)+tmp.row, 1]= ids[,2]
			idsAll[(1:tmp.row)+tmp.row, 2]= ids[,1]
			
			new.order = order(idsAll[,1])
			ids = NULL
			ids = idsAll[new.order,]
			
			uniqueIdx1 = unique(ids[,1])
			if(ifD) {
				print("unique idx:")
				print(uniqueIdx1)
			}
			beginIdx = match(uniqueIdx1, table=ids[,1])
			
			## rely on the natural order of ids[,2]
			id2 = ids[,2]
			ids = NULL
			
			senCt = length(uniqueIdx1)
			idProb = matrix(NA, nrow=senCt, ncol=3)
			idProb[,1]=uniqueIdx1
			idProb[,3]=beginIdx
			
			tmp = matrix(NA, nrow=senCt, ncol=2)
			tmp[,1]=hapProb[ uniqueIdx1 ]
			
			
			## to avoid processing too many if, leave the last one out for now
			
			tmp[1:(senCt-1), 2] = sapply(1:(senCt-1), FUN=function(i, idxcut, allmatch, prob){
						
						prob.idxR = (idxcut[i]):(idxcut[i+1]-1)
						prob.sum = sum(prob[  allmatch[ prob.idxR ]  ])
						return(prob.sum)
					}, idxcut=idProb[,3], allmatch=id2, prob=hapProb)
			
			## append the last one
			tmp[senCt, 2]= sum(hapProb[ id2[  (idProb[senCt,3]:length(id2))  ] ])
			
			idProb[,2]= tmp[,1]*tmp[,2]
			
			idProb[  (idProb[,2] > ((-1)*10^(-15))) &  (idProb[,2] < 0), 2]=0
			
			if(max(idProb[,2]<0)==1) stop("Sampling probability for strata A is wrong.")
			
			grpA = list(type="A", idProb=idProb, id2=id2, ct=senCt)
		}
		
		## also get D
		idProb=NULL
		len = length(hapProb)
		
		idProb = matrix(NA, nrow = len, ncol=2)
		idProb[,1] = 1:(len)
		
		tmp2 = rep(0, times=len)
		
		## if the first id is found in group A, we need to take out the probabilities
		if (grpA$ct > 0)  tmp2[ grpA$idProb[,1] ] = tmp[,2]
		
		idProb[,2] = hapProb*(1-hapProb-tmp2)
		tmp=NULL
		tmp2=NULL
		
		idProb[  (idProb[,2] > ((-1)*10^(-15))) &  (idProb[,2] < 0), 2]=0
		if(max(idProb[,2]<0)==1) stop("Sampling probability for strata D is wrong.")
		
		grpD = list(type="D", idProb=idProb, ct = len)
		
		return(list(grpA=grpA, grpD=grpD))
	}
	
	if (type == "B" | type=="C"){
		if(length(ids)==0) return(list(type=type, ct = 0))
		
		idProb = matrix(NA, nrow=length(ids) , ncol=2)
		# for type B or C, the ids are one column for homogenuous pair
		idProb[,1] = ids
		idProb[,2] = hapProb[ids]^2
		
		idProb[  (idProb[,2] > ((-1)*10^(-15))) &  (idProb[,2] < 0), 2]=0
		
		if(max(idProb[,2]<0)==1) stop("Sampling probability for strata B/C is wrong.")
		
		
		return(list(type=type, idProb=idProb, ct=length(ids)))
	}
	
}


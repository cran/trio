ESp.impuParent.heterKid <-
function(hapPair, hapProb, test=F){
	ifD = F
	fN = "ESp.impuParent.heterKid"
	hapCt = length(hapProb)
	
	## matrix to store the set info
	# first 2 col are hap idx, third is the probability
	setB = matrix(NA, ncol=3, nrow=5)
	setB[2,c(1,2)]=hapPair
	setB[3,c(1,2)]=rep(hapPair[1], 2)
	setB[5,c(1,2)]=rep(hapPair[2], 2)
	
	setB[2,3] = 2*hapProb[hapPair[1]]*hapProb[hapPair[2]] 
	setB[3,3] = (hapProb[hapPair[1]])^2
	setB[5,3] = (hapProb[hapPair[2]])^2
	setB[1,3] = 2*hapProb[hapPair[1]]*(1-sum(hapProb[hapPair]))
	setB[4,3] = 2*hapProb[hapPair[2]]*(1-sum(hapProb[hapPair]))
	
	
	sampleTb = matrix(NA, nrow=9, ncol=5)
	## first 2 col ae idx matched with setB
	sampleTb[,1]=c(1,1,1,2,2,2,3,3,3)
	sampleTb[,2]=c(4,2,5,4,2,5,4,2,5)
	## mating ratio
	sampleTb[,3]=c(2,2,2,2,1,2,2,2,2)
	## kids ratio
	sampleTb[,4]=c(.25, .25, .5, .25, .5, .5, .5, .5, 1)
	## last col is the sample prob
	
	## get the sample prob
	sampleTb[,5]= (setB[ sampleTb[,1], 3]) * (setB[ sampleTb[,2], 3]) * (sampleTb[,3]) * (sampleTb[,4])
	if(test){
		return(sampleTb)
	}
	if(ifD){
		print(paste(fN, "::sampling probability"))
		print(sampleTb)
	}
	
	sampleTb[,5]=sampleTb[,5]/sum(sampleTb[,5])
	chooseRow = sample(1:9, size=1, prob=sampleTb[,5])
	
	hap4idx = rep(NA, times=4)
	if(is.element(chooseRow, c(5,6,8,9))){
		#no need for resampling, the hapPair should come in as ordered
		hap4idx = c( setB[ sampleTb[chooseRow,1], c(1,2)], setB[ sampleTb[chooseRow,2], c(1,2)] )
		return(hap4idx)
	}else{
		hapProb[ hapPair ]=0
		sample.idx =  sample(1:hapCt, size=2, replace=T, prob=hapProb)
		# need to resampling A:A2
		if (chooseRow==1){
			hap4idx[ c(1,3) ] = hapPair
			hap4idx[ c(2,4) ] = sample.idx
			hap4idx = c(range(hap4idx[c(1,2)]), range(hap4idx[c(3,4)]))
			return(hap4idx)
		}
		# A:B
		if (chooseRow==2){
			hap4idx[ 1 ] = hapPair[1]
			hap4idx[ 2 ] = sample.idx[1]
			hap4idx[ c(3,4) ] = setB[2, c(1,2)]
			hap4idx = c(range(hap4idx[c(1,2)]), hap4idx[c(3,4)])
			return(hap4idx)
		}
		# A:C2
		if (chooseRow==3){
			hap4idx[ 1 ] = hapPair[1]
			hap4idx[ 2 ] = sample.idx[1]
			hap4idx[ c(3,4) ] = setB[5, c(1,2)]
			hap4idx = c(range(hap4idx[c(1,2)]), hap4idx[c(3,4)])
			return(hap4idx)
		}
		# B:A2
		if (chooseRow==4){
			hap4idx[ 3 ] = hapPair[2]
			hap4idx[ 4 ] = sample.idx[1]
			hap4idx[ c(1,2) ] = setB[2, c(1,2)]
			hap4idx = c(hap4idx[c(1,2)], range(hap4idx[c(3,4)]))
			return(hap4idx)
		}
		# C:A2
		if (chooseRow==7){
			hap4idx[ 3 ] = hapPair[2]
			hap4idx[ 4 ] = sample.idx[1]
			hap4idx[ c(1,2) ] = setB[3, c(1,2)]
			hap4idx = c(hap4idx[c(1,2)], range(hap4idx[c(3,4)]))
			return(hap4idx)
		}   
		
		
	}
	
	
}


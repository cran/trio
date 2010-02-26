ESp.impuParent.homoKid <-
function(hapPair, hapProb, test=F){
	
	hapCt = length(hapProb)
	
	## matrix to store the set info, no need anymore
	# first 2 col are hap idx, third is the probability
	
	sampleTb = matrix(NA, nrow=3, ncol=5)
	## first 2 col ae idx matched with setB
	sampleTb[,1]=c(1, 1, 2)
	sampleTb[,2]=c(1, 2, 2)
	## mating ratio, #not used anymore
	
	## kids ratio
	sampleTb[,4]=c(.25, .5, 1)
	## last col is the sample prob
	
	## get the sample prob
	
	## current approach: get the mating prob
	tp = hapProb[hapPair]
	matingProb = c(4*tp^2*(1-tp)^2, 4*tp^3*(1-tp), tp^4  )
	
	sampleTb[,5]= matingProb*sampleTb[,4]
	if(test){
		return(sampleTb)
	}
	sampleTb[,5]=sampleTb[,5]/sum(sampleTb[,5])
	chooseRow = sample(1:3, size=1, prob=sampleTb[,5])
	
	hap4idx = rep(NA, times=4)
	if(chooseRow==3){
		#no need for resampling
		hap4idx = rep(hapPair, times=4)
		return(hap4idx)
	}else{
		
		hapProb[ hapPair ]=0
		sample.idx =  sample(1:hapCt, size=2, replace=T, prob=hapProb)
		# need to resampling A:A2
		if (chooseRow==1){
			hap4idx[ c(1,3) ] = rep(hapPair, 2)
			hap4idx[ c(2,4) ] = sample.idx
			hap4idx = c(range(hap4idx[c(1,2)]), range(hap4idx[c(3,4)]))
			return(hap4idx)
		}
		# A:B
		if (chooseRow==2){
			hap4idx[ c(1,2,3) ] = rep(hapPair, 3)
			hap4idx[ 4 ] = sample.idx[1]
			hap4idx = c(hap4idx[c(1,2)],  range(hap4idx[c(3,4)]))
			return(hap4idx)
		}
	}
	
}


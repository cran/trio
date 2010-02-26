augBkFrame <-
function(hapBkMap, key, probLeftover = .01){
	keyIndex = which(hapBkMap$keys==key)
	bk = hapBkMap$bks[[keyIndex]]
	
	## extract the original exp and prob, then restandard prob
	estBkExp = as.character(bk[,hapBkMap$expCol])
	estBkProbRe = bk[,hapBkMap$probCol]*(1-probLeftover)
	
	## create the exhaust haplotypes 
	bkSnpLen = bk[1,hapBkMap$hapLenCol]
	exhaustExp = exhaustHapExp(lociCt=bkSnpLen, snpCoding=c(1,2))$hapStr
	exhaustExpCt = length(exhaustExp)
	
	## augmate the additional expression
	rematch = match(estBkExp, exhaustExp)
	
	resiProb = probLeftover / (exhaustExpCt-length(estBkExp))
	resiProb = rep(resiProb, times=exhaustExpCt)
	resiProb[rematch]=estBkProbRe
	
	base = bk[1,]
	newBk = NULL
	for( i in 1:exhaustExpCt ){
		newBk = rbind(newBk, base)
	}
	## replace the expression and probability
	newBk[,hapBkMap$probCol]=resiProb
	newBk[,hapBkMap$expCol]=as.integer(exhaustExp)
	
	## replace the necessary part in hapBkMap
	hapBkMap$bks[[keyIndex]]=newBk
	
	hapBkMap$bkLens[keyIndex]=exhaustExpCt
	
	## other parts, like df, dfStr,  in hapBkMap is not updated because it is not necessary
	
	return(hapBkMap)
	
}


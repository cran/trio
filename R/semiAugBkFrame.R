semiAugBkFrame <-
function(hapBkMap, key, probLeftOver = .01){
	
	ifD = F
	
	keyIndex = which(hapBkMap$keys==key)
	bk = hapBkMap$bks[[keyIndex]]
	
	## extract the original exp and prob, then restandard prob
	estBkExp = as.character(bk[,hapBkMap$expCol])
	estBkProbRe = bk[,hapBkMap$probCol]/(sum( bk[,hapBkMap$probCol] ))*(1-probLeftOver)
	
	## create the exhaust haplotypes 
	bkSnpLen = bk[1,hapBkMap$hapLenCol]
	exhaustExp = exhaustHapExp(lociCt=bkSnpLen, snpCoding=c(1,2))$hapStr
	exhaustExpCt = length(exhaustExp)
	
	## match the expression in the short list to the exhaustive list
	rematch = match(estBkExp, exhaustExp)
	
	## replace the expression's probability
	bk[,hapBkMap$probCol]=estBkProbRe
	bk$augIdx = rematch
	bk$resiProb = rep(probLeftOver, times=nrow(bk))
	
	hapBkMap$bks[[keyIndex]] = bk
	
	hapBkMap$augIdxCol = ncol(bk)-1
	hapBkMap$resiProbCol = ncol(bk)
	
	#if(ifD) print(bkFrame)
	
	## other parts, like df, dfStr,  in hapBkMap is not updated because it is not necessary
	return(hapBkMap)
	
}


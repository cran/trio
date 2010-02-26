dfToHapBkMap <-
function(data, keyCol=NULL, chCol, blockCol, expCol, probCol, hapLenCol, beginCol=NULL, endCol=NULL, snpBase=1, re.bf = T, re.javaGUI = T){
	
	if(is.null(keyCol)){
		data = cbind(data, key=paste(data[,chCol], data[,blockCol], sep="-"))
		keyCol = ncol(data)
	}
	
	allKeys = unique(data[,keyCol])
	keyCt = length(allKeys)
	
	bks = NULL
	snpLen = 0
	bkLens = NULL
	bkSnpLens = NULL
	markers_b = NULL
	markers_e = NULL
	
	df = data[,c( keyCol, chCol, blockCol, expCol, probCol, hapLenCol, beginCol, endCol)]
	
	if(!is.null(beginCol)){
		for( i in 1:keyCt){
			r = df[data[,keyCol]==allKeys[i], ]
			bks = c(bks, list(r))
			bkLens = c(bkLens, dim(r)[1])
			bkSnpLens = c(bkSnpLens, r[1, 6])
			snpLen = snpLen + r[1,6]
			markers_b = c(markers_b, r[1, 7])
			markers_e = c(markers_e, r[1, 8])
		}
	}else{
		for( i in 1:keyCt){
			r = df[data[,keyCol]==allKeys[i], ]
			bks = c(bks, list(r))
			bkLens = c(bkLens, dim(r)[1])
			bkSnpLens = c(bkSnpLens, r[1, 6])
			snpLen = snpLen + r[1, 6]
		}
	}
	
	
	dfStr = cbind(as.character(data[,keyCol]),
			as.character(data[,expCol]),
			as.character(data[,probCol]))
	
	dfMaster = cbind(as.character(allKeys),
			as.character(bkLens),
			as.character(bkSnpLens),
			
			as.character(markers_b),
			as.character(markers_e))
	
	hapBkMap = NULL
	if(re.bf){
		hapBkMap = c(hapBkMap, list(df=df))
	}
	
	if(re.javaGUI){
		hapBkMap = c(hapBkMap, list(dfStr=dfStr,       dfStrDim    = c(nrow(dfStr), ncol(dfStr)),
						dfMaster=dfMaster, dfMasterDim = c(nrow(dfMaster), ncol(dfMaster)),
						bkCumIndex=cumsum(bkLens)))
	}
	
	hapBkMap = c(hapBkMap, list(
					bks = bks,
					bkLens = bkLens,
					keys = as.character(allKeys),
					bkSnpLens = bkSnpLens, 
					snpBase =snpBase, snpLen = snpLen, keyCol=1, chCol=2, blockCol=3,
					expCol = 4, probCol = 5, hapLenCol=6))
	
	if(!is.null(beginCol)){
		hapBkMap = c(hapBkMap, list(markers_b = markers_b, markers_e = markers_e,
						beginCol=7, endCol=8))
	}
	
	return(hapBkMap)
}


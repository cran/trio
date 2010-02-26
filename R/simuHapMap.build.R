simuHapMap.build <-
function(hapInfoFrame=NULL, genoInfoFrame){
	
	## make it can take .csv file
	if(!is.null(hapInfoFrame)){
		if(is.character(hapInfoFrame)){
			hapInfoFrame = read.csv(hapInfoFrame, header=T, sep=",", as.is=T)
		}
	}
	
	## make it can take .csv file
	if(is.character(genoInfoFrame)){
		genoInfoFrame = read.csv(genoInfoFrame, header=T, sep=",", as.is=T)
	}
	
	## check input data
	test = match(c("prekey", "seq", "freq"), colnames(genoInfoFrame))
	if (sum(is.na(test))>=1)
		stop("Input argument, genoInfoFrame, does not have the required columns, i.e., prekey, seq, and freq.") 
	
#  if(is.null(hapInfoFrame)){
#    warning("NULL value is provided for input argument, hapInfoFrame.")
#  }
	
	genoMap.info = genoInfoFrame[, match(c("prekey", "seq", "freq"), colnames(genoInfoFrame))]
	
	genoMap = dfToGenoMap(df = genoMap.info, dataHeaders = c("prekey", 
					"seq", "freq"), genotype = c("11", "12", "22"), snpBase = 1)
	
	
	if(is.null(hapInfoFrame)){
		#print("here")
		
		onlyGeno = genoInfoFrame
		## build the df first
		key = paste( onlyGeno[,1], onlyGeno[,2], sep="-")
		key = unique(key)
		
		allGenoBk = NULL
		for( i in 1:(length(key))){
			x.b = (i-1)*3+1
			x.e = i*3
			oneBk = geno2hapFreq(prekey=onlyGeno[x.b,1], block=onlyGeno[x.b,2],
					genoFreq = onlyGeno[x.b:x.e, 3], snpBase=1,
					snpSeq = i)
			allGenoBk = rbind(allGenoBk, oneBk[,c(1,4,5)])
		}
		colnames(allGenoBk)=c("key","haploytype","frequency")
		#print(allGenoBk)
		
		simuMap = bkMap.constr(data=allGenoBk, keyCol=1, hapLenCol=NULL, expCol=2, probCol=3)
		simuSetup = list(NULL)
		simuSetup = c(simuSetup, newDf=list(allGenoBk), simuMap=list(simuMap)) 
	}else{
		bkFrame = hapInfoFrame
		bkFrame$key = paste(bkFrame[, 1], bkFrame[, 2], sep = "-h")
		
		hapBkMap = dfToHapBkMap(data = bkFrame, keyCol = ncol(bkFrame), 
				chCol = 1, blockCol = 2, expCol = 3, probCol = 4, 
				hapLenCol = 5, beginCol = 6, endCol = 7, snpBase = 1, 
				re.bf = T, re.javaGUI = T)
		hapBkGenoMap = bindHapBkGenoMaps(hapBkMap = hapBkMap, genoMap = genoMap)
		simuSetup = hapBkGenoMap2HapMap(hapBkGenoMap=hapBkGenoMap, reSimuMap=T)
	}
	return(simuSetup)
}


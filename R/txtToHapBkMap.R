txtToHapBkMap <-
function(txtF, delim=",", dataHeaders=NULL, sorted=F, ...){
	##txtF = "tblBlock.csv"
	##dataHeaders = NULL
	
	if(is.null(dataHeaders)){
		## assume that txtF has the first row as column names
		data = read.csv(file=txtF, header = T, sep=delim, na="missing")
	}else{
		## assume that column names of the data is passed from outside
		data = read.csv(file=txtF, header = F, sep=delim, na="missing")
		colnames(data) = dataHeaders
	}
	
	m = match(c("ch", "block", "hap", "freq", "hapLen","markers_b", "markers_e"), 
			colnames(data), 0)
	
	## assume except for markers_b and marekers_e, other variable should be presented
	if(min(m[1:5])==0) stop("One or more required variable(s) missing")
	
	## cast the data into the right format
	for ( i in m ){
		if(class(data[,i])=="factor") data[,i]=as.numeric(as.character(data[,i]))
	}
	
	## create key as a combination of chromosome and block
	key = paste(data[,m[1]], data[,m[2]], sep="-")
	
	
	df = cbind(key, data[,m])
	
	if(!sorted){
		df = df[ order(df$ch, df$block),]
	}
	
	hapBkMap = NULL
	if(m[6]==0){
		hapBkMap = dfToHapBkMap(df, keyCol=1, chCol=2, blockCol=3,
				expCol=4, probCol=5, hapLenCol=6, ...)
	}else{
		hapBkMap = dfToHapBkMap(df, keyCol=1, chCol=2, blockCol=3,
				expCol=4, probCol=5, hapLenCol=6,
				beginCol=7, endCol=8, ...)
	}
	return(hapBkMap)             
}


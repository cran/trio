txtToGenoMap <-
function(txtF, delim=",", dataHeaders=NULL, genotype=c("11", "12", "22"), snpBase=1){
	##txtF = "tblBlock.csv"
	##dataHeaders = NULL
	
	if(is.null(dataHeaders)){
		## assume that txtF has the first row as column names
		data = read.csv(file=txtF, header = T, sep=delim, as.is=T)
	}else{
		## assume that column names of the data is passed from outside
		data = read.csv(file=txtF, header = F, sep=delim, as.is=T)
		colnames(data) = dataHeaders
	}
	
	genoMap = dfToGenoMap(df=data, dataHeaders=dataHeaders,  genotype=genotype, snpBase=snpBase)
	
	return(genoMap)
	
}


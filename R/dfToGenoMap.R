dfToGenoMap <-
function(df, dataHeaders=NULL, genotype=c("11", "12", "22"), snpBase=1){
	##txtF = "tblBlock.csv"
	##dataHeaders = NULL
	
	if(is.null(dataHeaders)){
		## assume that txtF has the first row as column names
		data = df
	}else{
		## assume that column names of the data is passed from outside
		data = df
		colnames(data) = dataHeaders
	}
	
	m = match(c("prekey", "seq", "freq", "genotype"), colnames(data), 0)
	
	## assume except for markers_b and marekers_e, other variable should be presented
	if(min(m[2:3])==0) stop("One or more required variable(s) missing")
	
	
	if( !is.element("prekey", colnames(data)) ){
		## assume the seq of homo, hetero, homo
		data$prekey = I(rep("prekey", times=nrow(data)))
	}
	if( !is.element("genotype", colnames(data)) ){
		## assume the seq of homo, hetero, homo
		data$genotype = I(rep(genotype, times=nrow(data)/3))
	}
	mheader=c("prekey", "seq", "freq", "genotype")
	
	m = match(c("prekey", "seq", "freq", "genotype"), colnames(data), 0)
	## cast the data into the right format
	for ( i in 1:4 ){
		if(!is.na(m[i])){
			if (i==1){
				if(class(data[,m[i]])=="factor") data[,i]=as.character(data[,m[i]])
			}
			if (i==4){
				if(class(data[,m[i]])=="factor") data[,i]=as.numeric(as.character(data[,m[i]]))
			}
			if( i==2 | i==3){
				if(class(data[,m[i]])=="factor") data[,i]=as.numeric(as.character(data[,m[i]]))
			}
		}
	}
	
	chLastSNPIndex = util.sql.groupby(data[,m], groupCols=1,  varCol=2, type=c("max"))
	
	genoMap = list(df = data[,m], chLastSNPIndex=chLastSNPIndex, chCol=1, genomeSeqCol=2, probCol=3, expCol=4, snpBase=snpBase)
	return(genoMap)
	
}


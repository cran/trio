freq.build <-
function(hap=NULL, geno){
	
	## make it can take .csv file
	if(!is.null(hap)){
		if(is.character(hap)){
			hap = read.csv(hap, header=T, sep=",", as.is=T)
		}
	}
	
	## make it can take .csv file
	if(is.character(geno)){
		geno = read.csv(geno, header=T, sep=",", as.is=T)
	}
	
	## check input data
	test = match(c("prekey", "seq", "freq"), colnames(geno))
	if (sum(is.na(test))>=1)
		stop("Input argument, geno, does not have the required columns, i.e., prekey, seq, and freq.") 
	
	if(is.null(hap)){
		warning("NULL value is provided for input argument, hap.")
	}
	
	freq = c(genoMap.info = list(geno[, match(c("prekey", "seq", "freq"), colnames(geno))]),
			hapMap.info = list(hap))
	
	return(freq)
}


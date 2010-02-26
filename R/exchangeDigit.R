exchangeDigit <-
function(ma, cols=NULL, dig1Code=c(0, 1, 3, 2), dig2Code = c(0, 1, 2), action=c("1to2", "2to1")){
	ifD = F
	
	if(is.null(cols)){
		cols = c(1, ncol(ma))
	}
	
	fun.error = F
	if(length(cols)!=2) fun.error=T
	if(cols[2]>ncol(ma)) fun.error=T
	
	outputColCt = cols[2]-cols[1]
	if( outputColCt <=0) fun.error=T
	
	if(length(action)>=2) {
		
		if(length(action)>1) warning(paste("Two or more actions is request (", paste(action, collapse="; "),
							"). Only the first requested is performed."), sep="")
		action=action[1]
		
		if ((action!="1to2") & (action!="2to1"))
			stop(paste("Requested action, (",  action, "), is not implemented.", sep="" ))
	}
	ma.change = ma[, cols[1]:cols[2], drop=F]
	## internally, do not use NA to represent missing,
	code.0 = 0
	code.012 = dig2Code
	code.0123 = dig1Code
	
	if( max( is.na(dig2Code[2:3]))==1) stop("NA is not allow to represent the non-missing allele!")
	if( max( is.na(dig1Code[2:4]))==1) stop("NA is not allow to represent the non-missing genotype!")
	
	if(action=="1to2"){
		if( is.na(dig1Code[1]) ){
			code.123 = dig1Code[2:4]
			## if the min is the same as the default code for NA(=0), then use the min -1 
			if( max(code.123==0)==1) code.0 = min(code.123)-1
			code.0123 = c(code.0, code.123)
			ma.change = apply(ma.change, 1:2, FUN= util.vec.replace, orignal = dig1Code, replaceBy=code.0123)
		}
		outputColCt = outputColCt/2
		tmpBridge = dig2Code[c(1, 2, 2, 3)]
		tmpBridge2 = dig2Code[c(1, 2, 3, 3)]
		if(ifD) print(paste("tmpBridge:", paste(tmpBridge, collapse=";")))
		if(ifD) print(paste("tmpBridge2:", paste(tmpBridge2, collapse=";")))   
	}
	if(action=="2to1"){
		if( is.na(dig2Code[1]) ){
			#need to replace the given NA with 0 or others
			code.12 = dig2Code[2:3]
			if( max(code.12==0)==1) code.0 = min(code.12)-1
			code.012 = c(code.0, code.12)
		}
		ma.change = apply(ma.change, 1:2, FUN= util.vec.replace, orignal = dig2Code, replaceBy=code.012)
		outputColCt = outputColCt*2
		sumBase = matrix(c(2, 0, 0,
						0, 2, 0,
						0, 1, 1,
						0, 0, 2), byrow=F, nrow=3)
		sumCode = as.vector(matrix(code.012, ncol=3)%*%sumBase)
		if(ifD) print(paste("sumCode:", paste(sumCode, collapse=";")))
	}
	
	
	if (fun.error){
		stop("Argu, cols, refer to the ranges of the index for the columns in argu, ma.
						The current values are not valid.")
	}
	
	
	
	if(cols[1]>1) {
		ma.kept1 = ma[, 1:(cols[1]-1)]
	}else{
		ma.kept1 = NULL
	}
	
	ma.kept2=NULL
	if(cols[2]<ncol(ma)) ma.kept2 = ma[, (cols[2]+1):ncol(ma)]
	
	ma.ex=NULL
	if(action=="1to2"){
		
		
		ma.2dig1 = apply(ma.change, 1:2, FUN= util.vec.replace, orignal = code.0123, replaceBy=tmpBridge)
		ma.2dig2 = apply(ma.change, 1:2, FUN= util.vec.replace, orignal = code.0123, replaceBy=tmpBridge2)
		
		## make sure the small digit is at the front
		if(dig2Code[2]<dig2Code[3]){
			ma.ex = util.matrix.col.shuffle2(ma.2dig1, ma.2dig2)
		}else{
			ma.ex = util.matrix.col.shuffle2(ma.2dig2, ma.2dig1)
		}
		
	}
	
	if(action=="2to1"){
		
		seq1 = seq.int(from=1, to=ncol(ma.change), by=2)
		
		ma.sum = matrix(ma.change[, seq1]+ma.change[,seq1+1], ncol=length(seq1), byrow=F)
		
		## convert to Qing's 1-digit coding system. 3 is for hetero
		
		ma.ex = apply(ma.sum, 1:2, FUN= util.vec.replace, orignal = sumCode, replaceBy=dig1Code)
		
	}
	if(!is.null( ma.kept1)){
		all = cbind(ma.kept1, ma.ex)
	}else{
		all = ma.ex
	}
	if(!is.null( ma.kept2)){
		all = cbind(all, ma.kept2)
	}
	return(all)
	
}


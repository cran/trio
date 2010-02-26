bkMap.constr <-
function(data, keyCol, hapLenCol=NULL, expCol, probCol, alleleCode=1:2, ...){
	fN = "bkMap.constr"
	
	if( is.character(data)){
		data = read.csv(data, ...)
		#print(qp(fN, ": Info on the ", data, " file."))
		#print(str(data))
	}
	
	colNum = ncol(data)
	#print(keyVal)
	keyVal = unique(data[,keyCol]) # not reordered, following the original order
	#print(keyVal)
	
	bks = NULL
	snpLen = 0
	bkLens = NULL
	keys = NULL
	bkEndingIdx = NULL
	if(is.null(hapLenCol)){
		# if the input dataset does NOT contain info. on haplotype block length
		for( i in 1:length(keyVal)){
			if(i==1){
				r = data[data[,keyCol]==keyVal[i], ]
				
				## normalize the prob and check the haplen
				tmp.expLen  = unlist(lapply(r[,expCol], FUN=function(i){nchar(as.character(i))}  ))
				snpCt1 = tmp.expLen[1]
				if (min(tmp.expLen==snpCt1)==0) stop("At least one of the block has different lengths for haplotypes.")
				
				
				tmp.prob = r[,probCol]
				if (max(tmp.prob<0)==1) stop("At least one of the block has negative haplotype frequencies.")
				if (sum(tmp.prob)==0) stop("At least one of the block has haplotype frequencies sum up to 0.")          
				tmp.prob = tmp.prob/sum(tmp.prob)          
				r[,probCol]=tmp.prob
				
				bkLens = dim(r)[1]
				snpCtVec = rep(snpCt1, bkLens)
				
				r = cbind(r, snpCtVec)
				
				bks = list(r)          
				snpLen = snpLen + snpCt1
				keys = as.character(r[1,keyCol])
				snpCt = snpCt1
			}else{
				
				r = data[data[,keyCol]==keyVal[i], ]
				
				## normalize the prob and check the haplen
				tmp.expLen  = unlist(lapply(r[,expCol], FUN=function(i){nchar(as.character(i))}  ))
				snpCt1 = tmp.expLen[1]
				if (min(tmp.expLen==snpCt1)==0) stop("At least one of the block has different lengths for haplotypes.")
				
				tmp.prob = r[,probCol]
				if (max(tmp.prob<0)==1) stop("At least one of the block has negative haplotype frequencies.")
				if (sum(tmp.prob)==0) stop("At least one of the block has haplotype frequencies sum up to 0.")
				tmp.prob = tmp.prob/sum(tmp.prob)
				r[,probCol]=tmp.prob
				
				bkLens = c(bkLens, dim(r)[1])
				snpCtVec = rep(snpCt1, dim(r)[1])
				
				r = cbind(r, snpCtVec)
				
				bks = c(bks, list(r))
				snpLen = snpLen + snpCt1
				keys = c(keys, as.character(r[1,keyCol]))
				snpCt = c(snpCt, snpCt1)
			}
			
			hapLenCol = dim(r)[2]
		}
	}else{
		# if the input dataset contain info. on haplotype block length
		for( i in 1:length(keyVal)){
			if(i==1){
				r = data[data[,keyCol]==keyVal[i], ]
				
				## normalize the prob and check the haplen
				tmp.expLen  = unlist(lapply(r[,expCol], FUN=function(i){nchar(as.character(i))}  ))
				snpCt1 = tmp.expLen[1]
				if (min(tmp.expLen==snpCt1)==0) stop("At least one of the block has different lengths for haplotypes.")
				if (min(r[,hapLenCol]==snpCt1)==0) stop("At least one of the block provides inconsistant haplotype size.")
				
				tmp.prob = r[,probCol]
				if (max(tmp.prob<0)==1) stop("At least one of the block has negative haplotype frequencies.")
				if (sum(tmp.prob)==0) stop("At least one of the block has haplotype frequencies sum up to 0.")
				tmp.prob = tmp.prob/sum(tmp.prob)
				
				r[,probCol]=tmp.prob
				
				bks = list(r)
				bkLens = dim(r)[1]
				snpLen = snpLen + snpCt1
				keys = as.character(r[1,keyCol])
				snpCt = snpCt1
			}else{
				r = data[data[,keyCol]==keyVal[i], ]
				
				## normalize the prob and check the haplen
				tmp.expLen  = unlist(lapply(r[,expCol], FUN=function(i){nchar(as.character(i))}  ))
				snpCt1 = tmp.expLen[1]
				if (min(tmp.expLen==snpCt1)==0) stop("At least one of the block has different lengths for haplotypes.")
				if (min(r[,hapLenCol]==snpCt1)==0) stop("At least one of the block provides inconsistant haplotype size.")
				
				tmp.prob = r[,probCol]
				if (max(tmp.prob<0)==1) stop("At least one of the block has negative haplotype frequencies.")
				if (sum(tmp.prob)==0) stop("At least one of the block has haplotype frequencies sum up to 0.")
				tmp.prob = tmp.prob/sum(tmp.prob)
				
				r[,probCol]=tmp.prob
				
				bks = c(bks, list(r))
				bkLens = c(bkLens, dim(r)[1])
				snpLen = snpLen + snpCt1
				keys = c(keys, as.character(r[1,keyCol]))
				snpCt = c(snpCt, snpCt1)
			}
		}
	}
	bkMap = list(bks = bks, bkLens = bkLens, keys = keys, snpLen = snpLen, expCol = expCol, hapLenCol=hapLenCol, probCol = probCol, snpCt=snpCt, alleleCode=alleleCode )
	return(bkMap)
	
}


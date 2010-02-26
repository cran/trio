linkageFile.proc <-
function(data, snpIdxRange=NULL, key.prefix="", bk.sizes=NULL, action = c("outputTrio", "formTrio1digit", "missingReport", "Mendelian check", "freq estimate"), dig2Code=0:2, dig1Code=c(0,1,3,2), ... ){
	
	#dig2Code=0:2
	#dig1Code=c(NA,0,1,2)
	
	pedCol =1
	memCol=2
	affectCol = 6
	dadCol=3
	momCol=4
	txt.affect=2
	
	sep=""
	header=F
	
	
	trioTwoDigit = snpPREFileMatchTrio(txtF=data, sep=sep, header=header,
			pedCol =pedCol, memCol=memCol, affectCol =affectCol,
			dadCol=dadCol, momCol=momCol, txt.affect=txt.affect,
			logF=NULL)
	#print(dim(trioTwoDigit))
	
	if(is.null(snpIdxRange)) snpIdxRange = c(7, ncol(trioTwoDigit))
	
	snpStartLeftIndex=snpIdxRange[1]
	snpEndRightIndex =ncol(trioTwoDigit) - snpIdxRange[2] + 1
	
	re = list()
	if(is.element("outputTrio", action)){
		re = c(re, trio2digit=list(trioTwoDigit))
	}
	
	genos  = trioTwoDigit[, snpStartLeftIndex:(ncol(trioTwoDigit)-snpEndRightIndex+1)]
	snpNum = ncol(genos)/2
	
	
	## check linkage file code
	linkqing = unique(as.vector(unlist(genos)))
	tt.check = match(linkqing, dig2Code) 
	if(max(is.na(tt.check))==1)
		stop(paste("dig2Code is wrong.", "Existing codes in data:[",
						paste(linkqing, collapse=";", sep=""), "].")
		)
	
	snp1digit = exchangeDigit(ma=genos,
			cols=c(1,snpNum*2), dig1Code=dig1Code, dig2Code =dig2Code, action=c("2to1"))
	
	snp1digit.inside = exchangeDigit(ma=genos,
			cols=c(1,snpNum*2), dig1Code=c(0, 1, 3, 2), dig2Code =dig2Code, action=c("2to1"))
	
	trio1digit = cbind( trioTwoDigit[, c(pedCol, memCol)], snp1digit)
	if(is.element("formTrio1digit", action)){
		
		re = c(re, trio=list(trio1digit))
	}
	
	if(is.element("missingReport", action)){
		## check necessary parameters
		if( is.element(snpStartLeftIndex, c(pedCol, memCol, affectCol, dadCol, momCol))){
			warning("Cannot report missing information. The argument, snpIdxRange, includes one of the special column for linkage file.")
		}else{
			missSNPPos = findMissing(df=trioTwoDigit, is.1digit=F, snpStartLeftIndex=snpStartLeftIndex,
					snpEndRightIndex=snpEndRightIndex, dig1Code=NULL, dig2Code=dig2Code )
			re = c(re, missIdx = list(missSNPPos))
		}
	}
	
	if(is.element("Mendelian check", action)){
		snpTrio = matrix(snp1digit.inside, nrow=3, byrow=F)
		
		MedErr = matrix(NA, ncol=4, nrow=ncol(snpTrio))
		colnames(MedErr)=c("y", "x", "trio", "SNP")
		MedErr.ct = 0
		trioCt = nrow(snp1digit)/3
		
		## CHANGED: 2009: report the index in the trio1digit
		tmpDigit = 1
		#print(str(snpTrio))
		for ( i in 1:ncol(snpTrio)){
			tryCatch({
						tt = checkMendelianError(codedSNPTrio=snpTrio[,i], snpCoding=c(0,1,2,3))
					}, error = function(e){
						#print(paste("trioCt=", trioCt))
						#print(paste("snpStartLeftIndex=", snpStartLeftIndex))
						MedErr.ct <<- MedErr.ct +1
						#print(paste("MedErr.ct=", MedErr.ct, " i=", i))
						tttx = ceiling(i/trioCt)
						ttty = i%%trioCt
						if(ttty==0) ttty = trioCt
						## CHANGED: 2009: report the index in the trio1digit
						MedErr[MedErr.ct,] <<- c( (ttty-1)*3+1,  (tttx-1)*tmpDigit+3,  ttty,    tttx)
						#print( MedErr[MedErr.ct,,drop=F] )
					})          
		}
		
		if(MedErr.ct==0) {
			MedErr=NULL
			#print("No Mendelian error.")
			re = c(trio=list(trio1digit), MedErr=list(MedErr[1:MedErr.ct,,drop=F]))
		}else{
			#print("Found Mendelian error(s).")
			re = c(MedErr=list(MedErr[1:MedErr.ct,,drop=F]), trio.err=list(trio1digit))
		}
	}
	
	if(is.element("freq estimate", action)){
		## check necessary parameters
		if( is.element(snpStartLeftIndex, c(pedCol, memCol, affectCol, dadCol, momCol))){
			warning("Cannot provide frequencies estimation. The argument,  snpIdxRange, includes one of the special column for linkage file.")
		}else{
			
			tmp = ncol(trioTwoDigit)
			parTwoDigit = getBackParentGeno(trioDf=trioTwoDigit, famCol=pedCol, memCol=memCol,
					snpIdx=snpStartLeftIndex:(tmp-snpEndRightIndex+1), re.child=F, prefix=NULL)
			#print("###")
			#print(dim(parTwoDigit))
			
			if( is.null(bk.sizes) ){
				warning("Cannot provide frequencies estimation. The argument, bk.sizes, is not provided for user-specify option.")
			}else{
				
				map=freqmap.reconstruct(data=parTwoDigit, cols=c(3, ncol(parTwoDigit)), loci.ct=bk.sizes, is.1digit=F,
						dig1Code=NULL, dig2Code = dig2Code, key.prefix=key.prefix, start.base=1, ...)
				
				re = c(re, freq=list(map))
			}
			
		} ##if( is.element(snpStartLeftIndex, c(pedCol, memCol, affectCol, dadCol, momCol))){
		
	}
	
	return(re)
	
}


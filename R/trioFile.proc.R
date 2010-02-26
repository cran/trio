trioFile.proc <-
function(data, key.prefix="", bk.sizes=NULL, dig2Code=0:2, dig1Code=c(0,1,3,2),
		action = c("missingReport", "Mendelian check", "freq estimate"), ... ){
	
	
	sep=""
	header =F
	pedCol=1
	memCol=2
	snpIdxRange = c(3, ncol(data))
	is.1digit=T
	
#   if(!is.null(txtF)){
#     if(is.character(txtF)){
#          ## by default the text file has no header 
#          #trioTwoDigit = read.csv(file=data, header = header, sep=sep)
#          stop("The argument, data, cannot be a string.")
#     }else{
#          trioTwoDigit = data
#     }
#   }
	
	trioTwoDigit = data
	
	snpStartLeftIndex=snpIdxRange[1]
	snpEndRightIndex =ncol(trioTwoDigit) - snpIdxRange[2] + 1
	
	
	if(!is.1digit){
		genos  = trioTwoDigit[, snpStartLeftIndex:(ncol(trioTwoDigit)-snpEndRightIndex+1)]
		## check input, if is.1digit=T, dig2Code must be c(0,1,2), if is.1digit=F, dig1Code must be c(NA, 0,1,2,)
		code.Factor = apply( genos, 2, FUN=is.factor)
		if(max(code.Factor)==1) stop("Genotypes cannot be factors.")
		
		code.Factor = apply( genos, 2, FUN=is.numeric)
		if(min(code.Factor)==0) stop("Genotypes must be numerical.")
		
		code.Check = sort(unique(unlist(genos)))
		
		if( min(code.Check==dig2Code)==0 ){
			stop("Trio data is in 2-digit coding, but alleles are not represented by 0, 1, and 2.")
		}
		
		
		snpNum = ncol(genos)/2
		snp1digit = exchangeDigit(ma=genos,
				cols=c(1,snpNum*2), dig1Code=dig1Code, dig2Code =dig2Code, action=c("2to1"))
		
	}else{
		genos  = trioTwoDigit[, snpStartLeftIndex:(ncol(trioTwoDigit)-snpEndRightIndex+1)]
		
		## check input, if is.1digit=T, dig2Code must be c(0,1,2), if is.1digit=F, dig1Code must be c(NA, 0,1,2,)
		code.Factor = apply( genos, 2, FUN=is.factor)
		if(max(code.Factor)==1) stop("Genotypes cannot be factors.")
		
		code.Factor = apply( genos, 2, FUN=is.numeric)
		if(min(code.Factor)==0) stop("Genotypes must be numerical.")
		
		## check the coding for trio
		code.Check = unique(unlist(genos))
		
		code.Left = code.Check[!is.na(code.Check)]
		
		code.Match = match(code.Left, dig1Code[2:4])
		if( sum(!is.na(code.Match)) != length(code.Match)){
			stop("Non-missing genotypes are not coded as 0, 1, and 2")
		}
		
		snpNum = ncol(genos)
		snp1digit = genos
		## change to inside Qing's coding scheme, now assume the input code is c(NA, 0,1,2), no change, and no checking
		
#      if(min(dig1Code==c(0,1,2,3))==0){
#        snp1digit = apply(snp1digit, 1:2,  FUN= util.vec.replace, orignal = dig1Code, replaceBy=c(0,1,2,3))
		
	}
	
	snp1digit.inside = apply(snp1digit, 1:2,  FUN= util.vec.replace, orignal = dig1Code,
			replaceBy=c(0,1,3,2))
	
	re = list()
	
#   if(is.element("formTrio1digit", action)){
#     trio1digit = cbind( trioTwoDigit[, c(pedCol, memCol)], snp1digit)
#     re = c(re, trio1digit=list(trio1digit))
#   }
	
	if(is.element("missingReport", action)){
		## check necessary parameters
		if( is.element(snpStartLeftIndex, c(pedCol, memCol))){
			warning("Cannot report missing information. The argument, snpStartLeftIndex, is one of the special column for the trio file.")
		}else{
			
			if(is.1digit){
				#print(here)
				#print(snp1digit.inside)
				#print(snpStartLeftIndex)
				#print(snpEndRightIddex)
				missSNPPos = findMissing(df=snp1digit.inside, is.1digit=T, snpStartLeftIndex=snpStartLeftIndex-2,
						snpEndRightIndex=snpEndRightIndex, dig1Code=c(0,1,3,2), dig2Code=dig2Code)
			}else{
				missSNPPos = findMissing(df=trioTwoDigit, is.1digit=F, snpStartLeftIndex=snpStartLeftIndex,
						snpEndRightIndex=snpEndRightIndex, dig1Code=NULL, dig2Code=dig2Code )
			}
			
			if(is.null(missSNPPos)){
				
				re = c(re, missIdx = list(NULL))
			}else{
				re = c(re, missIdx = list(missSNPPos))
			}
		}
	}
	
	if(is.element("Mendelian check", action)){
		tmpDigit=2
		if(is.1digit) tmpDigit = 1
		snpTrio = matrix(snp1digit.inside, nrow=3, byrow=F)
		
		MedErr = matrix(NA, ncol=4, nrow=ncol(snpTrio))
		colnames(MedErr)=c("y", "x", "trio", "SNP")
		MedErr.ct = 0
		trioCt = nrow(snp1digit)/3
		
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
						MedErr[MedErr.ct,] <<- c( (ttty-1)*3+1,  (tttx-1)*tmpDigit+1+snpStartLeftIndex-1,  ttty,    tttx)
						#print( MedErr[MedErr.ct,,drop=F] )
					})          
		}
		
		if(MedErr.ct==0) {
			MedErr=NULL
			#print("No Mendelian error.")
		}else{
			#print("Found Mendelian error(s).")
		}
		re = c(re, MedErr=list(MedErr[1:MedErr.ct,,drop=F]))
		
	}
	
	if(is.element("freq estimate", action)){
		## check necessary parameters
		if(is.null(bk.sizes)) warning("Cannot provide frequencies estimation. The argument, bk.sizes, is missing.")
		
		if( is.element(snpStartLeftIndex, c(pedCol, memCol))){
			warning("Cannot provide frequencies estimation. The argument, snpStartLeftIndex, is one of the special column for linkage file.")
		}else{
			
			
			tmp = ncol(trioTwoDigit)
			parTwoDigit = getBackParentGeno(trioDf=trioTwoDigit, famCol=pedCol, memCol=memCol,
					snpIdx=snpStartLeftIndex:(tmp-snpEndRightIndex+1), re.child=F, prefix=NULL)
			#print(dim(parTwoDigit))
			
			if( is.null(bk.sizes) ){
				warning("Cannot provide frequencies estimation. The argument, bk.sizes, is not provided for user-specify option.")
			}else{
				
				map=freqmap.reconstruct(data=parTwoDigit, cols=c(3, ncol(parTwoDigit)), loci.ct=bk.sizes, is.1digit=is.1digit,
						dig1Code=dig1Code, dig2Code = dig2Code, key.prefix=key.prefix, start.base=1, ...)
				
				re = c(re, freq=list(map))
			}
			
		} ##if( is.element(snpStartLeftIndex, c(pedCol, memCol, affectCol, dadCol, momCol))){
		
	}
	
	return(re)
	
}


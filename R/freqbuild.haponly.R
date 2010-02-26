freqbuild.haponly <-
function(hap, alleleCode = 1:2){
	## only hap is provided.
	## need to find out the hap with 1 loci, and get the geno freq from the 1 loci
	ifD = F
	
	## check the maximum length for imputation
	bkMap = bkMap.constr(data=hap, keyCol=1, hapLenCol=NULL, expCol=2, probCol=3, alleleCode=alleleCode )
	if(sum(bkMap$snpCt >=8 )>0)  stop("At least one of the block size exceeds the maximum value of 7 loci.") 
	
	prekey = rep("ch", nrow(hap))
	cumidx = cumsum(bkMap$snpCt)
	markers.e = cumidx
	
	markers.b = c(1, (cumidx+1)[-length(cumidx)])
	others = lapply(  1:(length(bkMap$keys)), FUN=function(i, rep1, rep2, rep3, reptimes){
				outData1 = rep(rep1[i], times=reptimes[i])
				outData2 = rep(rep2[i], times=reptimes[i])
				outData3 = rep(rep3[i], times=reptimes[i])
				ooo = cbind(outData1, outData2, outData3)
			},
			rep1=bkMap$snpCt, rep2=markers.b, rep3=markers.e, reptimes=bkMap$bkLens)
	
	others.ma =NULL
	for( dd in others){
		others.ma = rbind(others.ma, dd)
	}
	others = matrix(others, nrow=nrow(hap), ncol=3, byrow=T)
	
	hapFrame = cbind(prekey, hap, others.ma)
	colnames(hapFrame)=c("prekey", "block", "hap", "freq", "hapLen", "markers_b", "markers_e")
	
	## need to remove the singletones
	hapFrame = hapFrame[ (hapFrame[,5]!=1), ]
	
	
	genoFrame = data.frame(prekey=rep("ch", bkMap$snpLen*3),
			seq=rep(1:bkMap$snpLen, each=3),
			freq=rep(NA, bkMap$snpLen*3))
	
	single = which(bkMap$snpCt == 1)
	if(length(single)>0){
		
		## find singletons.
		for( ss in single ){
			freq.info = bkMap$bks[[ss]]
			allele.freq=freq.info[,3][match(freq.info[,2], alleleCode)]
			geno1 = allele.freq[1]^2
			geno2 = 2*allele.freq[1]*allele.freq[2]
			geno3 = allele.freq[2]^2
			genoFrame[ ((ss-1)*3+1):(ss*3) , 3]=c(geno1, geno2, geno3)
		}
		
		
	}
	if(nrow(hapFrame)>0){
		freq.haponly = freq.build(hap=hapFrame, geno=genoFrame)
	}else{
		freq.haponly = freq.build(hap=NULL, geno=genoFrame)
	}
	if(ifD) {
		print(hapFrame)
		print(dim(genoFrame))
		print(genoFrame[1:6,])
	}
	return(freq.haponly)
}


freqmap.reconstruct <-
function(data, cols=NULL, loci.ct, is.1digit=F, dig1Code=c(0, 1, 2, 3), dig2Code = c(0, 1, 2), key.prefix="", start.base=1, ...){
   ifD = F

   ## is NA is error used to represent the missing, change it to 0, or the max(allele code)+1
   
   if (is.1digit){
     if(is.null(cols)) cols=c(1, ncol(data))
     data = exchangeDigit(ma=data, cols=cols, dig1Code=dig1Code,
       dig2Code = dig2Code, action=c("1to2")) [,(cols[1]: ((cols[2]-cols[1]+1)*2+cols[1]-1))]
     miss.code = dig1Code[1]
     #str(data)
   }else{
     data = data[, (cols[1]:cols[2])]
     miss.code = dig2Code[1]
   }
   
   ## need to construct the haplotype freq file
   ## with hap key, haplotype, haplotype prob 
   ct = length(loci.ct)
   geno.code = as.vector(matrix(dig2Code[2:3], ncol=2) %*% matrix(c(2,0, 1,1, 0,2), ncol=3, byrow=F))

   #if(sum(loci.ct>7)>=1) stop("Cannot process haplotype block with 8 or more loci.")
   if(sum(loci.ct)!=(ncol(data)/2)) {
     if(ifD) print(sum(loci.ct))
     if(ifD) print(ncol(data)/2)
     stop("Dismatching number of loci with the number of column in data.")

   }

   hapMap.info= NULL
   genoMap.info=NULL
   
   startId = start.base
   
   ## process the file if only one block exist
   if(ct==1){

     if(loci.ct>=2){
       hap.re = haplo.em(geno=data, ...)
       
       hap.prob = hap.re$hap.prob
    
       hap.exp =  hap.re$haplotype
    
       hap.expVec = as.integer(apply(hap.exp, 1, paste, collapse=""))
       #m = match(c("ch", "block", "hap", "freq", "hapLen","markers_b", "markers_e"), 
       df = data.frame( prekey= rep(key.prefix, length(hap.prob)),
                   block =rep(1, length(hap.prob)),
                   hap=as.vector(hap.expVec),
                   freq = hap.prob,
                   hapLen = rep(loci.ct, length(hap.prob)),
                   markers_b=rep(startId, length(hap.prob)),
                   markers_e=rep(startId+loci.ct-1, length(hap.prob))
                  )
       hapMap.info=df
     }

     dd = data[,1]+data[,2]
       
     tb = table(factor(unlist(dd), levels=geno.code))
     tb.freq = tb/sum(tb)

     #ch	seq	freq	genotype
     df = data.frame(prekey= I(rep(key.prefix, 3)),
                   seq =rep(startId, 3),
                   freq = as.vector(tb.freq)
                  )
     genoMap.info=df
     
     
     return(re=list(hapMap.info=hapMap.info, genoMap.info=genoMap.info))
     
   }

    ## process the file if only more than one block exist
   idx.end = cumsum(loci.ct)*2
   idx.start = c(1, (idx.end+1)[-length(loci.ct)])
   cut.rg = cbind(idx.start, idx.end)


   all.df.key2 = NULL
   all.df.oth2 = NULL

   snp.b = startId
   snp.e = startId-1
   tmp1 = 1
   
   allsnp = T
   for( i in 1:ct){
      keys =NULL
      oth=NULL
      snp.b = snp.e+1
      snp.e = snp.b+loci.ct[i]-1
      if(loci.ct[i]!=1){
		  allsnp = F
          ## for Hap
          hap.re = haplo.em(geno=data[, ((cut.rg[i,1]):(cut.rg[i,2]))], ...)

          # print(hap.re)
          hap.prob = hap.re$hap.prob
       
          hap.exp  = hap.re$haplotype
       
          hap.expVec = as.integer(apply(hap.exp, 1, paste, collapse=""))

          keys = rep(key.prefix, length(hap.prob))
          
          oth =cbind(block =  rep(i, length(hap.prob)),
                   hap=hap.expVec,
                   freq = hap.prob,
                   hapLen = rep(loci.ct[i], length(hap.prob)),
                   markers_b=rep(snp.b, length(hap.prob)),
                   markers_e=rep(snp.e,length(hap.prob))
                  )

          all.df.key2 = c(all.df.key2, keys)
          all.df.oth2 = rbind(prekey=all.df.oth2, oth)
          #print(all.df.oth2[1:3,])
          
      } #if(loci.ct[i]!=1){
    } #for( i in 1:ct){

	if(!allsnp){
	    rownames(all.df.key2)=NULL
	    rownames(all.df.oth2)=NULL
	    hapMap.info = data.frame(prekey=I(all.df.key2), all.df.oth2)
	}else{
		hapMap.info = NULL
	}
 
    ## gen one genomap for every SNP
    dd.t = data[, (min(cut.rg)):(max(cut.rg))]
    # convert into two columns
    dd.ma = matrix(unlist(dd.t), ncol=2, byrow=F)
    dd = dd.ma[,1]+dd.ma[,2]

    ## remove NA first
    #print(table(dd))
    ## convert into n columns each column is for one SNP
    dd = matrix(dd, nrow=nrow(data), byrow=F)
    tb = apply(dd, 2, FUN= function(col, geno.code){
      tt = table(factor(unlist(col), levels=geno.code))
      tt.freq = tt/sum(tt)
      tt.freq
    }, geno.code=geno.code)

    #print(tb[1:3,])

    #ch	seq	freq	genotype
 
    genoMap.info = data.frame(prekey= I(rep(key.prefix,  3*ncol(dd))),
        seq =rep(startId : (startId + cumsum(loci.ct)[length(loci.ct)] -1), each=3),
        freq = as.vector(tb) )

   	return(re=list(hapMap.info=hapMap.info, genoMap.info=genoMap.info))
   
}


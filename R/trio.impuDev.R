trio.impuDev <-
function(trio1digit, freqMaps, bk.seq=NULL, dig1Code=0:3, subInF = "PPC", dir=NULL, job=1, prefix="", impu.missingOnly=T){

  if( subInF!="PPC" ) stop("Subject in the data should follow Father, Mother and child sequence.")
  
  if(is.character(trio1digit) ) {
    data = read.csv(trio1digit, header=F)
  }else{
    data = trio1digit
  }
  
  toolboxNames = toolbox.load(freqMaps=freqMaps)
  
  if(!is.null(  get(toolboxNames$freqMap)$hapBkOnlyMap )){
	  bk.ssize =  get(toolboxNames$freqMap)$hapBkOnlyMap$bkSnpLens
	  if (max(bk.ssize>=8)==1) stop("At least one haplotype block has 8 or more SNPs in the block. Method fails.")	  
  }  
  
  
  # inside the function, assume 0 1 3 2 for NA, homo, hetero, homo, and 0 1 2 for NA, allele1, allele2
  dig1Default = c(0, 1, 3, 2)
  ## if the given code for 1-digit is not (0 1, 3, 2), change it.

  if( sum(dig1Default == dig1Code)!=4 ){
    data.geno = apply(data[, c(-1, -2)], 1:2, FUN= util.vec.replace, orignal = dig1Code, replaceBy=dig1Default)
    data = cbind(data[, 1:2], data.geno)
    
  }
  ##TODO!!! confirm the number
  #print(toolboxNames)
  
  bkCt = nrow(get(toolboxNames$freqMap)$genomeMarkerInfo)

  #print(data[1:3, 1:10])
  if(is.null(bk.seq)) bk.seq = 1:bkCt
  if(impu.missingOnly){

    ## get the hapPair without saving.
    imputed = impuBk.scheduler(raw=data, idx=bk.seq, job=job,
                toolname=toolboxNames,
                freqMaps=NULL, dir=dir,
                is.1digit=T, dig1Code=dig1Default, dig2Code=0:2,
                reType=F,
                reHap=qp(prefix, "impuHap_bk", bk.seq[1]),
                logF=NULL,
                logErr=qp(prefix, "impuErr_bk", bk.seq[1])
      )
    if( sum(dig1Default == dig1Code)!=4 ){
      imputed = apply(imputed, 1:2, FUN= util.vec.replace, orignal = dig1Default, replaceBy=dig1Code)
    }
     return(imputed)
  }

  ## get the hapPair without saving.
  imputed = impuBkTDT.scheduler(raw=data, idx=bk.seq, job=job,
                toolname=toolboxNames,
                freqMaps=NULL, dir=dir,
                is.1digit=T, dig1Code=dig1Default, dig2Code=0:2,
                reType=F,
                reHap=qp(prefix, "impuHap_bk", bk.seq[1]),
                logF=NULL,
                logErr=qp(prefix, "impuErr_bk", bk.seq[1])
  )

  if( sum(dig1Default == dig1Code)!=4 ){
     imputed = apply(imputed, 1:2, FUN= util.vec.replace, orignal = dig1Default, replaceBy=dig1Code)
  }
  return(imputed)
 
}


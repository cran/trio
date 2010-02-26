trio.impu <-
function(triodd, freq,  impu.missingOnly=T){
  dir=NULL

  dig1Code=c(NA,0,1,2)


  ## check number of loci match with map or not
  if( (ncol(triodd)-2) != nrow(freq$genoMap.info)/3 )
    stop("The number of SNPs in the frequency file does not equal the number of SNPs in trio dataset.")

  
  #if( subInF!="PPC" ) stop("Subject in the data should follow Father, Mother and child sequence.")
  
  if(is.character(triodd) ) {
    data = read.csv(triodd, header=F)
  }else{
    data = triodd
  }
  code.Factor = apply( data[,c(-1,-2)], 2, FUN=is.factor)
  if(max(code.Factor)==1) stop("Genotypes cannot be factors.")

  code.Factor = apply( data[,c(-1,-2)], 2, FUN=is.numeric)
  if(min(code.Factor)==0) stop("Genotypes must be numerical.")
  
  ## check the coding for trio
  code.Check = unique(unlist(data[,c(-1,-2)]))
  
  if( sum(is.na(code.Check))==0 ){
    if (impu.missingOnly) {
      print("No missing genotype. Return original data.")
      return(triodd[,-c(1,2)])
    }
  }

  code.Left = code.Check[!is.na(code.Check)]

  code.Match = match(code.Left, dig1Code[2:4])
  if( sum(!is.na(code.Match)) != length(code.Match)){
    stop("Genotypes are not coded as 0, 1, nor 2")
  }

  toolboxNames = toolbox.load(freqMaps=freq)
  
  if(!is.null(  get(toolboxNames$freqMap)$hapBkOnlyMap )){
	  bk.ssize =  get(toolboxNames$freqMap)$hapBkOnlyMap$bkSnpLens
	  if (max(bk.ssize>=8)==1) stop("At least one haplotype block has 8 or more SNPs in the block. Method fails.")	  
  }

  # inside the function, assume 0 1 3 2 for NA, homo, hetero, homo, and 0 1 2 for NA, allele1, allele2
  dig1Default = c(0, 1, 3, 2)
  ## if the given code for 1-digit is not (0 1, 3, 2), change it.

  if (is.na(dig1Code[1])){
    ttt = data[,c(-1,-2)]
    ttt[is.na(ttt)]=max(dig1Code, na.rm=T)+1

    data[,c(-1,-2)]=ttt
    dig1Code[1]=max(dig1Code, na.rm=T)+1
  }

  
  if( sum(dig1Default == dig1Code)!=4 ){
    data.geno = apply(data[, c(-1, -2)], 1:2, FUN= util.vec.replace, orignal = dig1Code, replaceBy=dig1Default)
    data = cbind(data[, 1:2], data.geno)
  }
  ##TODO!!! confirm the number
  #print(toolboxNames)
  
  bkCt = nrow(get(toolboxNames$freqMap)$genomeMarkerInfo)

  if(impu.missingOnly){

    ## get the hapPair without saving.
    imputed = impuBk.scheduler(raw=data, idx=1:bkCt, job=1,
                toolname=toolboxNames,
                freqMaps=NULL, dir=dir,
                is.1digit=T, dig1Code=dig1Default, dig2Code=0:2,
                reType=F, reHap=NULL, logF=NULL, logErr="")
    if( sum(dig1Default == dig1Code)!=4 ){
      imputed = apply(imputed, 1:2, FUN= util.vec.replace, orignal = dig1Default, replaceBy=dig1Code)
    }
     return(imputed)
  }

  #print("get impuBkTDT.scheduler")
  ## get the hapPair without saving.
  imputed = impuBkTDT.scheduler(raw=data, idx=1:bkCt, job=1,
                toolname=toolboxNames,
                freqMaps=NULL, dir=dir,
                is.1digit=T, dig1Code=dig1Default, dig2Code=0:2,
                reType=F, reHap=NULL, logF=NULL, logErr="")

  if( sum(dig1Default == dig1Code)!=4 ){
     imputed = apply(imputed, 1:2, FUN= util.vec.replace, orignal = dig1Default, replaceBy=dig1Code)
  }
  return(imputed)
 
}


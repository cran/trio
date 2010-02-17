impuBk.scheduler <-
function(raw, idx, job=1, toolname=NULL, freqMaps=NULL, dir="", is.1digit=T, dig1Code, dig2Code, reType=F, reHap=NULL, logF=NULL, logErr=""){

  fStr ="[impuBk.scheduler:]"
  ifD = F

  if(!is.null(dir)){
    if(dir==""){
      print(paste("Imputation information files(s) will be saved under current working directory:",  getwd()))
    }else{
      print(paste("Create directory:", dir, ", where imputation information file(s) are saved.", sep=""))
      dir.create(path=dir, showWarnings = TRUE)
      if (!is.null(logF)){
        logF = file.path(dir, qp(logF, ".txt"))
      }
      if (!is.null(reHap)){
        reHap = file.path(dir, reHap)
      }
      logErr = file.path(dir, logErr)
    }
  }else{
    if(job==1) {
#      print(
#          paste("Value for argument dir is NULL. Function will not save any information about imputation. In case of error, log file will be saved under current working directory:", getwd()))
#       paste("Function will not save any information about imputation. In case of error, log file will be saved under current working directory:", getwd()))
      reHap = NULL
      logF = NULL
    }else{
      if (!is.null(logF)){
        logF = qp(logF, ".txt")
      }
#      print(
#          paste("Value for argument dir is NULL. However, users ask to do multiple imputation. ",
#                "Imputation information files(s) will be saved under current working directory:",  getwd()))
    }
  }

  # inside the function, assume 0 1 2 3 for NA, homo, homo, heter and 0 1 2 for NA, allele1, allele2

  ## first check the existance of required parameters
  if(is.null(toolname)){
    if(is.null(freqMaps)){
      stop("Values for both arguments, toolboxName and freqMaps, are missing.")
    }else{
      toolname = toolbox.load(freqMaps=freqMaps)
    }
  }

  ## find the start and end position for the blocks of interest

  all.genomeMarkerInfo = get(toolname$freqMap)$genomeMarkerInfo

  bd.start = all.genomeMarkerInfo[idx[1], 2, drop=T]
  bd.end   = all.genomeMarkerInfo[idx[ length(idx) ], 3, drop=T]

  snpOffset =bd.start-1

  ## change digit to 1 digit coding using Qing's coding scheme
  if(!is.1digit){
      ## excluding the parents cols
      genos  = raw[, 2+( ((bd.start-1)*2+1):(bd.end*2) ),drop=F  ]
      snpNum = ncol(genos)/2
      snp1digit = exchangeDigit(ma=genos,
        cols=c(1,snpNum*2), dig1Code=c(0, 1, 3, 2), dig2Code =dig2Code, action=c("2to1"))
  }else{
      ## excluding the parents cols
      genos  = raw[, 2+(bd.start:bd.end), drop=F]
      snpNum = ncol(genos)
      snp1digit = genos
      #dig1Code.inside = dig1Code[c(1, 2, 4, 3)]
      
      ## change to inside Qing's coding scheme
      if(min(dig1Code==c(0,1,3,2))==0){
        snp1digit = apply(snp1digit, 1:2,  FUN= util.vec.replace, orignal = dig1Code, replaceBy=c(0,1,3,2))
      }
  }

  trioCt = nrow(snp1digit)/3
  maxRow = trioCt*job

  ## 18 columns: bk, trio, x1, x2, y1, y2, hap_f, hap_m, hap_c1-4
  if(reType){
    imputBkRecord = matrix(NA, ncol = 19, nrow=maxRow)
    impuDummy = imputBkRecord 
  }else{
    imputBkRecord = matrix(NA, ncol = 18, nrow=maxRow)
    impuDummy = imputBkRecord 
  }
  imputBkRecord.ct = 0
  errorTrap = NULL
  
  tryCatchEnv = new.env(parent=baseenv())
  assign("trapID", 0, env=tryCatchEnv)
  assign("errorTrap", errorTrap, env=tryCatchEnv)

  if(!is.null(get(toolname$freqMap)$hapBkOnlyMap)){
	  all.hapIndex = get(toolname$freqMap)$hapIndex
	  
	  hapBkOnlyMap.vars=list()
	  tmp.hapMap =  get(toolname$freqMap)$hapBkOnlyMap
	  hapBkOnlyMap.vars$resiProbCol= tmp.hapMap$resiProbCol
	  hapBkOnlyMap.vars$augIdxCol= tmp.hapMap$augIdxCol
	  hapBkOnlyMap.vars$probCol= tmp.hapMap$probCol

  }else{
	  all.hapIndex = c(-1)
	  hapBkOnlyMap.vars=list()
  }
  snpCoding = 0:3
  snpBase = 0:2
  genoProb = genGenoProb()

    for( unit in idx){
      if (is.element(unit, all.hapIndex)){
        ## haplotype, for each block, search every trio for missingness
        ## find block boundary
        bk.bd = unlist(all.genomeMarkerInfo[unit, c(2,3)])
        bk.geno = snp1digit[, (bk.bd[1]:bk.bd[2])]

        snpCt = bk.bd[2]-bk.bd[1]+1

        exhaustHap = get(toolname$exp)
        exhaustHap = exhaustHap[1:2^snpCt, 1:snpCt]

        bk.genoRowComp = rowSums(bk.geno!=0) == snpCt 
        bk.genoTrioComp = matrix(bk.genoRowComp, ncol=3, byrow=T)
        bk.missTrioIdx = (1:trioCt) [ rowSums(bk.genoTrioComp) < 3 ]
    
        for (famId in bk.missTrioIdx){
          x1 = bk.bd[1]
          x2 = bk.bd[2]
          y1 = (famId-1)*3+1
          y2 = famId*3
          trioBlock = snp1digit[y1:y2, x1:x2 -  snpOffset]
          if(ifD) print( paste(fStr, " processing fam index:", famId))
          if(ifD) print(trioBlock)
          
          tryCatch({
            replace = ESp.imputBlock(appVarNames=toolname, 
                           trioBlock=trioBlock, snpLen=snpCt, bkIdx = unit, job=job,
                           snpCoding = snpCoding, snpBase=snpBase, reType = reType,  logF=logF, 
                           hapBkOnlyMap.vars=hapBkOnlyMap.vars
                            )

            imputBkRecord.ct = imputBkRecord.ct + 1
            imputBkRecord[((imputBkRecord.ct-1)*job+1):(imputBkRecord.ct*job) , 1:6 ] = matrix(
                           rep(c(unit, raw[y1, 1], y1, y2, x1, x2), times=job), ncol=6, byrow=T)
            
            imputBkRecord[((imputBkRecord.ct-1)*job+1):(imputBkRecord.ct*job) ,
                          7:ncol( imputBkRecord ) ] = replace
            
            fam.hapIdx1 = replace[1, c(1,3,5)]
            fam.hapIdx2 = replace[1, c(2,4,6)]
            hap.str1 = exhaustHap[fam.hapIdx1, ]
            hap.str2 = exhaustHap[fam.hapIdx2, ]

            ## convert string into digits
            geno.FMCMa = covDipStr2CodedGeno(hap.str1, hap.str2, subjectCt=3, snpLen=snpCt)
 
            if(ifD) print(geno.FMCMa)
            if(ifD) print(snp1digit[y1:y2, x1:x2-  snpOffset])

            if( sum(abs(geno.FMCMa - trioBlock)[ trioBlock!=0 ])!=0 ) stop("No matching")

            snp1digit[y1:y2, x1:x2-  snpOffset] = geno.FMCMa

            if(!is.null(logF)){
              logl(logF, paste("For trio #", y2/3, ", Choose hap idx=[",
                               paste(replace[1:6], collapse=".", sep=""),
                               "]"))
              logl(logF, paste("Choose hap exp=[",
                    paste( apply(geno.FMCMa, 1, FUN=paste, collapse=""), collapse=".", sep=""), "]",
                               sep=""))
              logl(logF, paste("missing bk data:",
                    paste( apply(trioBlock, 1, FUN=paste, collapse=""), collapse="."),
                               sep=""))
      
            }
            if(ifD){
              print("-------------------")
              print(  c(y1, y2, x1, x2)  )              
            }

#            if( (imputBkRecord.ct!=0) & (imputBkRecord.ct %%cutpt==0) & (!is.null(reHap)) ){  
#                  ## if there is record and upto certain number, need to be outputed
#                  ## evaluated after each trio
#                  write.table(imputBkRecord[1:(imputBkRecord.ct*job) ,,drop=F],
#                              file=paste(reHap, "_imputBkRecord.csv", sep=""), sep=",",
#                              append = T, row.names = F, col.names=F)
#                  imputBkRecord.ct = 0
#                  imputBkRecord = impuDummy
#             }                       
           },  error = function(e) {
               ##print(qTraceback())
               traceback()
			   b = get("trapID", env=tryCatchEnv)
			   b = b +1 
			   assign("trapID", b, env=tryCatchEnv)
               #trapID <<- trapID+1
	
               ## HARD CODE!!!HARD CODE: trio id is assumed to the be first one
			   #errorTrap <<- rbind(errorTrap, errorInfo)
	           errorInfo = c(trapID=b, bkIdx=unit, pedgree=raw[y1,1], case=raw[(y1+2),2], c(y1, y2, x1, x2))
			   a = get("errorTrap", env=tryCatchEnv)
			   a = rbind(a, errorInfo)
			   assign("errorTrap", a, env=tryCatchEnv)
			   
               
                if(!is.null(logF)){
                  logl(logF, paste("\nError trap id=(", b, ") and details for errors:", sep=""))
                  logl(logF, paste("Error block index {idx=", unit, "}------------", sep=""))
                  logl(logF, paste("Error trap famId=(", famId, ") and details for errors:", sep=""))
                  
                  logl(logF, paste("missing bk data:",
                               paste( apply(trioBlock, 1, FUN=paste, collapse=""), collapse="."),
                               sep=""))
                  
                  #logl(logF, errTrace)
                }         
               }, warn =function(w) {
                 print("ImpuBlock::Warnings")
                 
            } ) ## tryCatch

           if(!is.null(logF)){
             logl(logF, paste("end imputing block index {idx=", unit, "}------------\n", sep=""))
           }
              
        } ## for (famId in bk.missTrioIdx){
         if( (imputBkRecord.ct!=0) & (!is.null(reHap)) ){  
                  ## if there is record and upto certain number, need to be outputed
                  ## evaluated when the block is over
                 write.table(imputBkRecord[1:(imputBkRecord.ct*job) ,,drop=F],
                             file=paste(reHap, "_imputBkRecord.csv", sep=""), sep=",",
                             append = T, row.names = F, col.names=F)
         }
         imputBkRecord.ct = 0
         imputBkRecord = impuDummy
      }else{
        # print("A genotype")  ## genotype
        bk.bd = unlist(all.genomeMarkerInfo[unit, 2])
        tmp.Map =  get(toolname$freqMap)$genoOnlyMap
        popuProb =  tmp.Map$bks[[bk.bd]][,tmp.Map$probCol]

        snpCt = 1

        bk.geno = snp1digit[,bk.bd -  snpOffset]        
        bk.genoRowComp = as.integer(bk.geno!=0) 
        bk.genoTrioComp = matrix(bk.genoRowComp, ncol=3, byrow=T)
        bk.missTrioIdx = (1:trioCt) [ rowSums(bk.genoTrioComp) < 3 ]

        for (famId in bk.missTrioIdx){
           x1 = bk.bd
           x2 = bk.bd
           y1 = (famId-1)*3+1
           y2 = famId*3
           trioBlock = snp1digit[y1:y2, x1:x2 -  snpOffset]
           if(ifD) print( paste(fStr, " processing fam index:", famId))
           if(ifD) print(trioBlock)
           
           tryCatch({
 
             replace = imputGeno( trioBlock, job=job, genoProb=genoProb,
                    popuProb = popuProb, data.order="FMC", snpCoding=snpCoding)

             #print(replace)
             imputBkRecord.ct = imputBkRecord.ct + 1
             imputBkRecord[((imputBkRecord.ct-1)*job+1):(imputBkRecord.ct*job) , 1:6 ]  =  matrix(
                            rep(c(unit, raw[y1, 1], y1, y2, x1, x2), times=job), ncol=6, byrow=T)

             imputBkRecord[((imputBkRecord.ct-1)*job+1):(imputBkRecord.ct*job),
                           c(7, 9, 11, 13, 15, 17) ] = replace
             
             geno.FMCMa = replace[1, 1:3]
             if( sum(abs(geno.FMCMa - trioBlock)[ trioBlock!=0 ])!=0 ) stop("No matching")

             snp1digit[y1:y2, x1:x2 -  snpOffset] = replace[1,1:3]

#              if( (imputBkRecord.ct!=0) & (imputBkRecord.ct %%cutpt==0) & (!is.null(reHap)) ){  
#                   ## if there is record and upto certain number, need to be outputed
#                   ## evaluated after each trio                
#                 write.table(imputBkRecord[1:(imputBkRecord.ct*job) ,,drop=F],
#                             file=paste(reHap, "_imputBkRecord.csv", sep=""), sep=",",
#                             append = T, row.names = F, col.names=F)
#                 imputBkRecord.ct = 0
#                 imputBkRecord = impuDummy
#              }
             
           },  error = function(e) {
               ##print(qTraceback())
               traceback()
			   
			   b = get("trapID", env=tryCatchEnv)
			   b = b +1 
			   assign("trapID", b, env=tryCatchEnv)
			   #trapID <<- trapID+1
			   
	           ## HARD CODE!!!HARD CODE: trio id is assumed to the be first one
			   #errorTrap <<- rbind(errorTrap, errorInfo)
			   errorInfo = c(trapID=b, bkIdx=unit, pedgree=raw[y1,1], case=raw[(y1+2),2], c(y1, y2, x1, x2))
    		   a = get("errorTrap", env=tryCatchEnv)
			   a = rbind(a, errorInfo)
			   assign("errorTrap", a, env=tryCatchEnv)			   
     
                if(!is.null(logF)){
                  logl(logF, paste("\nError trap id=(", b, ") and details for errors:", sep=""))
                  logl(logF, paste("Error block index {idx=", unit, "}------------", sep=""))
                  logl(logF, paste("Error trap famId=(", famId, ") and details for errors:", sep=""))
                  
                  logl(logF, paste("missing bk data:",
                               paste(trioBlock, collapse="."),
                               sep=""))
                }
               }, warn =function(w) {
                 # print("ImpuBlock::Warnings")                 
            } ) ## tryCatch

           if(!is.null(logF)){
             logl(logF, paste("end imputing block index {idx=", unit, "}------------\n", sep=""))
           }
         } ## for (famId in bk.missTrioIdx){
         if( (imputBkRecord.ct!=0) & (!is.null(reHap)) ){  
                ## if there is record, write out after each singletonn SNP 
                write.table(imputBkRecord[1:(imputBkRecord.ct*job) ,,drop=F],
                            file=paste(reHap, "_imputBkRecord.csv", sep=""), sep=",",
                            append = T, row.names = F, col.names=F)
         }
         imputBkRecord.ct = 0
         imputBkRecord = impuDummy
     
      } ## if (is.element(unit, all.hapIndex)){
    } ## for( unit in idx){

	errorTrap=get("errorTrap", env=tryCatchEnv)
    if(!is.null(errorTrap)){
      write.table(errorTrap, file=paste(logErr, "_errorTrap.csv", sep=""), sep=",",
                 append = F, row.names = F, col.names = TRUE)
     
    }
    
    return (snp1digit)
}


trio.simu.proposed <-
function(bkMap, rule, caseNo, datasetCt=1, infoS="simuPropInfo", exInfoS="exSimuPropInfo", ddir=NULL, startIdx=NULL, spStrata.saveFN = NULL,  spStrata.name=NULL, reControl =F, dig1Code=0:3 ){

  if (is.null(ddir)) {
    infoS=NULL
    exInfoS = NULL
  }
  ifD = F
  
  info = bkMap.HRCB.LRCB.split(bkMap, rule)

  bkMapNS = info$bkMapNS

  varPrefix = paste(sample(letters[1:26], size=10), sep="", collapse="")
  varOriginalName = "stepstone"
  finalUse = paste(varPrefix, varOriginalName, sep="_")
  if(ifD) print(qp("finalUse=", finalUse))

  if(is.null(spStrata.saveFN)){
    ## if the spStrata is not saved earlier, need to generate
    if(!is.null(spStrata.name)){
      ## then we can choose to save it
      spStrata = bkMap.HRCB.Esp1Rule.Base(bkMap, rule, baseName=spStrata.name)
      if(ifD) print(qp("No stepstone object is previously saved. Generate and save the object as:", spStrata.name, ".RData"))

      #finalUse = spStrata.name
      assign(finalUse, spStrata, env=.GlobalEnv)
      
    }else{
      ## or not save it, just used. Not recommend.
      spStrata = bkMap.HRCB.Esp1Rule.Base(bkMap, rule, baseName=NULL)

      #finalUse = spStrata
      assign(finalUse, spStrata, env=.GlobalEnv)
      #print(qp("No stepstone object is previously saved. Generate but not save the object."))
    }
    
  }else{
    ## if the spStrata is saved earlier, just load the saved one
    
    if( is.character(spStrata.saveFN)){
      if(ifD) print(qp("Stepstone object is previously saved as:", spStrata.saveFN, ".RData. Do not need to generate it."))

      tryCatch({tmp = load(paste(spStrata.saveFN, ".RData", sep=""))
               assign(finalUse, get(tmp[1]), env=.GlobalEnv)  },
             error = function(e){
               stop(paste("Cannot open step-stone file '", spStrata.saveFN, ".RData'.", sep=""))
             })

    }else{
      #print(qp("Stepstone object is given."))
      #finalUse = spStrata.saveFN
      assign(finalUse, spStrata.saveFN, env=.GlobalEnv)
    }
  }
  
  preObj = bkMap.HRCB.Esp1Rule.genoSeq(bkMap, rule, re.probOnly = F)

  
  simuTrio = rep(list(NA), length=datasetCt)

  if(is.null(ddir) &  (datasetCt!=1))
    if(ifD) print( "Request to generate multiple datasets, and will return as a list.")

  for( i in 1:datasetCt){

        if(!is.null(infoS)){
          trioData1 = HRCB.Esp1Rule.spTrioOnBase(bkMap=NULL, preObj=preObj, spStrata=finalUse,
                                       rule=rule, caseNo,
                                       ifS = qp(infoS, "_", i) , reControl=reControl)
        }else{
          trioData1 = HRCB.Esp1Rule.spTrioOnBase(bkMap=NULL, preObj=preObj, spStrata=finalUse,
                                       rule=rule, caseNo,
                                       ifS = NULL , reControl=reControl)
        }
    #tryCatch({
        if(!is.null(exInfoS)){
          trioData2 = bkMap.LRCB.spTrio(bkMapNS, caseNo=caseNo, ifS = qp(exInfoS, "_", i), reControl=reControl)
        }else{
          trioData2 = bkMap.LRCB.spTrio(bkMapNS, caseNo=caseNo, ifS = NULL, reControl=reControl)
        }
     #}, warning=function(w){print(w)})
        
        simuTrioData = trioMerge(trioData1, trioData2, colName1=info$newColName, colName2=info$noCausalColName)

        ## get back to common coding scheme, 3 for heter
        if( is.null(dim(simuTrioData))){
          snp1recode =  util.vec.replace(simuTrioData, orignal = c(0,1,3,2), replaceBy= dig1Code)
        }else{
          snp1recode =  apply(simuTrioData, 1:2, FUN= util.vec.replace, orignal = c(0,1,3,2), replaceBy=dig1Code)
        }

        ## adding other stuff
        if(ifD) print("A")
        if(ifD) print(dim(snp1recode))
        famid = rep(1:caseNo, each=3)
        pid = rep(1:3, times=caseNo)
        trio.snp = cbind(famid, pid, snp1recode)
        colnames(trio.snp) = c("famid", "pid", paste("snp", 1:(ncol(snp1recode)), sep=""))
        
        if(is.null(ddir)){
          simuTrio[[i]]=trio.snp
        }else{

         if (ddir==""){
           if (is.null(startIdx)){
             save(trio.snp, file=qp("trioSimu_", i, ".RData"))
           }else{
             save(trio.snp, file=qp("trioSimu_", startIdx+i-1, ".RData"))
           }
         }else{
           if (is.null(startIdx)){
             save(trio.snp, file=file.path(ddir, qp("trioSimu_", i, ".RData")))
           }else{
             save(trio.snp, file=file.path(ddir, qp("trioSimu_", startIdx+i-1, ".RData")))
           }
         } 
       }

  }

  if(is.null(ddir)){
    return(simuTrio)
  }else{
    return(NULL)
  }

}


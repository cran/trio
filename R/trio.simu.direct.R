trio.simu.direct <-
function (bkMap, rule, caseNo, datasetCt=1, infoS="simuDirInfo", ddir=NULL, baseObj.saveFN=NULL, baseObj.name=NULL,  reControl=F, dig1Code=0:3 ){

  info = bkMap.HRCB.LRCB.split(bkMap, rule)

  bkMapNS = info$bkMapNS

  varPrefix = paste(sample(letters[1:26], size=10), sep="", collapse="")
  varOriginalName = "stepstone"
  finalUse = paste(varPrefix, varOriginalName, sep="_")
  #print(qp("finalUse=", finalUse))

  if(is.null(baseObj.saveFN)){
    ## if the spStrata is not saved earlier, need to generate
    if(!is.null(baseObj.name)){
      ## then we can choose to save it
     
      print(qp("No stepstone object is saved. Choose to generate and save the object as:", baseObj.name))
      
      matingTbInfo = bkMap.HRCB.famMap(info$bkMapS, rule, newColName=info$newColName,  ifS=infoS, baseName=baseObj.name)

      #finalUse =  baseObj.name
      assign(finalUse, matingTbInfo, env=.GlobalEnv)
      
    }else{
      ## or not save it, just used. Not recommend.
      print(qp("No stepstone object is saved. Choose to generate but not save the object."))
      matingTbInfo = bkMap.HRCB.famMap(info$bkMapS, rule, newColName=info$newColName,  ifS=infoS, baseName=NULL)

      #finalUse = matingTblInfo
      assign(finalUse, matingTbInfo, env=.GlobalEnv)
    }
    
  }else{
    ## if the spStrata is saved earlier, just load the saved one
   
    if( is.character(baseObj.saveFN)){
      print(qp("Stepstone object is saved before as:", baseObj.saveFN, ". Do not need to generate it."))
      tmp = load(paste(baseObj.saveFN, ".RData", sep=""))
      assign(finalUse, get(tmp[1]), env=.GlobalEnv)
      
    }else{
      print(qp("Stepstone object is given."))
      assign(finalUse, baseObj.saveFN, env=.GlobalEnv)
    }
  }

  
  simuTrio = rep(list(NA), length=datasetCt)

  if(is.null(ddir) &  (datasetCt!=1)){
    print( "Request to generate multipe datasets, but do not give path for save. Will return as a list.")
  }

  for( i in 1:datasetCt){
      if(is.null(infoS)){
        trioData1 = HRCB.famMap.spTrio(caseNo=caseNo, matingTbName=finalUse,
          ifS =NULL , reControl=reControl)
      }else{
        trioData1 = HRCB.famMap.spTrio(caseNo=caseNo, matingTbName=finalUse,
          ifS = qp(infoS, "_", i) , reControl=reControl)
      }
      
      trioData2 = bkMap.LRCB.spTrio(bkMapNS, caseNo=caseNo, reControl=reControl)
      
      simuTrioData = trioMerge(trioData1, trioData2, colName1=info$newColName, colName2=info$noCausalColName)

      if( is.null(dim(simuTrioData))){
        snp1recode =  util.vec.replace(simuTrioData, orignal = c(0,1,3,2), replaceBy= dig1Code)
      }else{
        snp1recode =  apply(simuTrioData, 1:2, FUN= util.vec.replace, orignal = c(0,1,3,2), replaceBy=dig1Code)
      }

      if(is.null(ddir)){
        simuTrio[[i]]=snp1recode
      }else{
         if (ddir==""){
           save(snp1recode, file=qp("trioDirSimu_", i, ".RData"))
         }else{
           save(snp1recode, file=file.path(ddir, qp("trioDirSimu_", i, ".RData")))
         } 
      }
    }

  if(is.null(ddir)){
    return(simuTrio)
  }else{
    return(NULL)
  }
}


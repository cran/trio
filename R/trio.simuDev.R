trio.simuDev <-
function(bkMap=NULL,  sigStr="g9=11 and g13=11", sigType, para=c(-1, .5),  caseNo=10, datasetCt=1,  dig1Code=0:3, ddF=NULL, startIdx=NULL, stepstone.saveFN=NULL, stepstone.name=NULL, spSupHap.namePrefix=NULL, verbose=T, reControl=F){

  print(paste("Try to simulated trio data: # datasets=", datasetCt, "; # trio =", caseNo, sep=""))


  if(!is.null(ddF)){
    if(ddF==""){
      print(paste("Simulated data file(s) and other result file(s) are saved under current working directory."))
      infoS.s = spSupHap.namePrefix
      exInfoS.s = paste(spSupHap.namePrefix, "LRCB", sep="")
    }else{
      print(paste("Create directory:", ddF, ", where simulated data file(s) and other result file(s) are saved.", sep=""))
      dir.create(path=ddF, showWarnings = TRUE)
      infoS.s = file.path(ddF, spSupHap.namePrefix)
      exInfoS.s = file.path(ddF, paste(spSupHap.namePrefix, "LRCB", sep=""))
    }

    if(is.null(spSupHap.namePrefix))  infoS.s = NULL
  }else{
    print(paste("Value for argument ddF is NULL. Function will return data as a list, ignoring input for argument, spSupHap.namePrefix"))
    infoS.s = NULL
    exInfoS.s = NULL
  }

  if(is.null(bkMap)){
    print(paste("Value for argument bkMap is NULL. Function will use package default object, simuBkMap."))
    data("simuBkMap")
    bkdata = get("simuBkMap")
    bkMap = bkMap.constr(data=bkdata, keyCol=1, hapLenCol=NULL, expCol=2, probCol=3, alleleCode=1:2)
  }

  #print("Info on data object for haplotype block frequencies:")
  #print(str(bkMap, max.level=1))

  ##  Dec08Change!!! allow D/R coding
  rule =signalRule.2strata.build(sigStr=sigStr, sigType=sigType, para=para)

  #print("Info on data object for risk factor, signalRule:")
  #print(rule$slist[[1]]$str)
  if(is.null(spSupHap.namePrefix)){
    infoS.s = NULL
    exInfoS.s = NULL
  }else{
    if(!is.null(ddF)){
      if(ddF==""){
        infoS.s = spSupHap.namePrefix
        exInfoS.s = paste(spSupHap.namePrefix, "LRCB", sep="")
      }else{
        infoS.s = file.path(ddF, spSupHap.namePrefix)
        exInfoS.s = file.path(ddF, paste(spSupHap.namePrefix, "LRCB", sep=""))
      }
    }else{
      infoS.s = NULL
      exInfoS.s = NULL
    }
  }

  ptm=proc.time()


  trioData = trio.simu.proposed(bkMap=bkMap,
                    rule=rule,
                    caseNo=caseNo,
                    datasetCt=datasetCt,
                    infoS=infoS.s,
                    exInfoS=exInfoS.s,
                    ddir=ddF,
                    startIdx=startIdx,
                    spStrata.saveFN = stepstone.saveFN,
                    spStrata.name = stepstone.name,
                    reControl =reControl,  dig1Code=dig1Code 
                    )
  
  if(verbose) print(paste("Time used to generate ", datasetCt,  " dataset(s).", sep=""))
  if(verbose) print(proc.time() - ptm)
  gc.e = gc()
 
  if(verbose) print("Info on memory usage:")
  if(verbose) print(gc.e)

  print(paste("Finished generating ", datasetCt,  " dataset(s).", sep=""))

  return(trioData)
}


trio.simuOLD <-
function(bkMap=NULL,  interaction="9R and 13R", alpha, beta, n=10, rep=1,  stepstone.saveFN=NULL, stepstone.name=NULL, verbose=T){

  sigType="D/R"
  para = c(alpha, beta)
  dig1Code=c(4,0,1,2)
  ddF = NULL
  spSupHap.namePrefix=NULL
  alleleCode=1:2

  print(paste("Try to simulated trio data: # datasets=", rep, "; # trio =", n, sep=""))
  
  if(!is.null(ddF)){
    if(ddF==""){
      print(paste("Simulated data file(s) and other result file(s) are saved under current working directory:",  getwd()))
      infoS.s = spSupHap.namePrefix
    }else{
      print(paste("Create directory:", ddF, ", where simulated data file(s) and other result file(s) are saved.", sep=""))
      dir.create(path=ddF, showWarnings = TRUE)
      infoS.s = file.path(ddF, spSupHap.namePrefix)
    }

    if(is.null(spSupHap.namePrefix))  infoS.s = NULL
  }else{
    #print(paste("Value for argument ddF is NULL. Function will return data as a list, ignoring input for argument, spSupHap.namePrefix"))
    infoS.s = NULL
  }

  if(is.null(bkMap)){
    print(paste("Value for argument bkMap is NULL. Function will use package default object, simuBkMap."))
    data("simuBkMap")
    bkdata = get("simuBkMap")
    bkMap = bkMap.constr(data=bkdata, keyCol=1, hapLenCol=NULL, expCol=2, probCol=3, alleleCode=alleleCode)
  }

  #print("Info on data object for haplotype block frequencies:")
  #print(str(bkMap, max.level=1))
   ##  Dec08Change!!! allow D/R coding
  rule =signalRule.2strata.build(sigStr=interaction, sigType=sigType, para=para)

  #print("Info on data object for risk factor, signalRule:")
  #print(rule$slist[[1]]$str)

  ptm=proc.time()

  trioData = trio.simu.proposed(bkMap=bkMap,
                    rule=rule,
                    caseNo=n,
                    datasetCt=rep,
                    infoS=infoS.s,
                    exInfoS=NULL,
                    ddir=ddF,
                    spStrata.saveFN = stepstone.saveFN,
                    spStrata.name = stepstone.name,
                    reControl =F,  dig1Code=dig1Code 
                    )
  
  if(verbose) print(paste("Time used to generate ", rep,  " dataset(s).", sep=""))
  if(verbose) print(proc.time() - ptm)
  gc.e = gc()
 
  if(verbose) print("Info on memory usage:")
  if(verbose) print(gc.e)

  print(paste("Finished generating ", rep,  " dataset(s).", sep=""))

  return(trioData)
}


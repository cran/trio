bkMap.HRCB.Esp1Rule.Base <-
function(bkMap, rule, baseName=NULL, dig1Code=0:3){
  ifD = F
  fn = "bkMap.HRCB.Esp1Rule.Base::"
  if(ifD) print(paste(fn, "begin"))

  if (length( rule$slist)>1) stop("This function work only with one SignalRule list")

  
  signalRule = rule$slist[[1]]

  HRCBIdx = bkMap.findHRCBIdx(bkMap, as.integer(signalRule$var.idx), re.keys=F, unique=T)
  superHRCBMap = bkMap.superHRCB(bkMap, uniBkIndexes=HRCBIdx, re.probOnly = T)
  
  # obtain super-hap and their probabilities
  supHapProb = superHRCBMap
  supLen = length(supHapProb)

    
  # if only one term in the rule, it doesn't matter whether the coef is positive or negative
  # if more than one term, need to choose the negative as the reference
  ## 20Dec08Change!!!
  HRCBGrps = bkMap.ESp.apply1Rule(bkMap, signalRule)

  if(ifD){
    print( HRCBGrps )
    #return(NULL)
  }

  ## construct the objects for four sampling stratum
  
  if (nrow(HRCBGrps$A)==0){
    ## if no heter pairs associated with high risk
    #print("check")
    #print(HRCBGrps$A)
  }
  
  HRCBStra.AD = HRCBSpGrp.cons(supHapProb, HRCBGrps$A, type="A")
  #print(HRCBStra.AD)

  HRCBStra.A = HRCBStra.AD$grpA

  HRCBStra.D = HRCBStra.AD$grpD

  #print(supHapProb)
  #print(HRCBGrps$B)
  HRCBStra.B = HRCBSpGrp.cons(supHapProb, HRCBGrps$B, type="B")
  #print(HRCBStra.B)
  
  if(ifD) {
    print(HRCBStra.A)
    print(HRCBStra.D)
    print(HRCBStra.B)
  }
  if(length( HRCBGrps$B )>0){
    tmp.Cidx = (1:supLen)[-HRCBGrps$B]
  }else{
    tmp.Cidx = (1:supLen)
  }
  HRCBStra.C = HRCBSpGrp.cons(supHapProb, tmp.Cidx, type="C")

  if(ifD){
    # check the sum of prob
    print(paste(fn, " check the sum of prob:"))
    print(sum(HRCBStra.A$idProb[,2],
        HRCBStra.B$idProb[,2],
        HRCBStra.C$idProb[,2],
        HRCBStra.D$idProb[,2]))
     print(HRCBStra.A)
     print(HRCBStra.B)
     print(HRCBStra.C)
     print(HRCBStra.D)
  }

  ## restandadize the prob
  straCt = c(HRCBStra.A$ct, HRCBStra.B$ct, HRCBStra.C$ct, HRCBStra.D$ct)
  straCumCt = cumsum(straCt)
  sttProb = c(HRCBStra.A$idProb[,2],
        HRCBStra.B$idProb[,2],
        HRCBStra.C$idProb[,2],
        HRCBStra.D$idProb[,2])
  sttProb = sttProb/sum(sttProb)

  #print(paste(fn, " risk allele freq = ", sum( sttProb [1: straCumCt[2] ])))

  
  if(HRCBStra.A$ct>0) HRCBStra.A$idProb[,2] = sttProb[1:straCumCt[1]]
  if(HRCBStra.B$ct>0) HRCBStra.B$idProb[,2] = sttProb[(straCumCt[1]+1): straCumCt[2] ]
  if(HRCBStra.C$ct>0) HRCBStra.C$idProb[,2] = sttProb[(straCumCt[2]+1): straCumCt[3] ]
  HRCBStra.D$idProb[,2] = sttProb[(straCumCt[3]+1): straCumCt[4] ]

  HRCBStra = list(HRCBStra.A = HRCBStra.A, HRCBStra.B = HRCBStra.B,
                  HRCBStra.C = HRCBStra.C, HRCBStra.D = HRCBStra.D)

  if(!is.null(baseName))  save(HRCBStra, file=paste(baseName, ".RData", sep=""))

  return(HRCBStra)

}


ESp.imputBlock <-
function(appVarNames,  trioBlock, snpLen=ncol(trioBlock),  bkIdx, job=1, snpCoding, snpBase, reType=F, logF=NULL,  hapBkOnlyMap.vars){

  ## TODO!!! making missed only disappear
  fStr ="[ESp.imputBlock:]"
  ifD = F
  # inside the function, assume snpCoding as c( 0, 1, 2, 3 )for NA, homo, homo, heter and snpBase as c(0, 1, 2) for NA, allele1, allele2

  
  if(ifD) print( paste(fStr, " processing block index:", bkIdx))

  if( min( c(snpCoding==c(0,1,2,3),  snpBase ==c(0,1,2))) <1 )  
     stop (paste("\nData configuration is not right:\n", "snpCoding=[", paste(snpCoding, collapse=";", sep=""),
                                                       "] snpBase=[", paste(snpBase, collapse=";", sep=""), "]", sep=""))
  allhapKeys = get(appVarNames$freqMap)$hapIndex
  semiMapFrame = get(appVarNames$freqMap)$hapBkOnlyMap$bks[[ match(bkIdx, allhapKeys)]]

  ## the third row is the child
  child = trioBlock[3,]
  father = trioBlock[1,]
  mother  = trioBlock[2,]

  compMissing.trio = as.logical(apply(trioBlock, 1, sum)==0)

  cHapBkInfoMap = hapGenoBlockProc(child,  snpCoding=snpCoding)
  reqIn =  cHapBkInfoMap$homoIn
  reqDig = cHapBkInfoMap$homoDigit
  

  ## if no parents is completely missing
  if( !(compMissing.trio[1] | compMissing.trio[2])   ){
      
      fHapBkInfoMap = hapGenoBlockProc(father,  snpCoding=snpCoding)
      mHapBkInfoMap = hapGenoBlockProc(mother,  snpCoding=snpCoding)
    
    
      fHapFiltered = procSemiAugMap(appVarNames, fHapBkInfoMap, snpLen)
      mHapFiltered = procSemiAugMap(appVarNames, mHapBkInfoMap, snpLen)
      
      ## obtain completed parent filtered dip, knock off the one doesn't match with homozygous index
      fDipTb = fHapBkIdx2DipTb(appVarNames, fHapBkInfoMap, idxList=fHapFiltered, snpCoding,
                                      reqIn=reqIn, reqDigits = reqDig, expression=fHapBkInfoMap$ori, snpLen)
     
      mDipTb = fHapBkIdx2DipTb(appVarNames, mHapBkInfoMap, idxList=mHapFiltered, snpCoding,
                                      reqIn=reqIn, reqDigits = reqDig, expression=mHapBkInfoMap$ori, snpLen)


    
    if (compMissing.trio[3]){

      if(ifD) print( "no parents is compMiss, child is compMiss, sample parents.")
      if(ifD) print( "Type 2: F(d)M(d)C(CM)")
      ## no parents is compMiss, child is compMiss, sample parents.

      ## sample the non missing parent
      fProb = exDipProbSemiAugMap(fDipTb, semiMapFrame=semiMapFrame,
                                  hapBkOnlyMap.vars$resiProbCol, hapBkOnlyMap.vars$augIdxCol, hapBkOnlyMap.vars$probCol, snpLen)
      ## sample the non missing parent
      mProb = exDipProbSemiAugMap(mDipTb, semiMapFrame=semiMapFrame,
                                  hapBkOnlyMap.vars$resiProbCol, hapBkOnlyMap.vars$augIdxCol, hapBkOnlyMap.vars$probCol, snpLen)

      trioHapIdx = matrix(NA, nrow=job, ncol=13)
      trioHapIdx[,13]=rep(2, job)
      for( ss in 1:job){
         ## sample fa then sample mo...
         chooseRow = sample(length(fProb), size=1, prob=fProb)
         tmpFather = fDipTb[chooseRow,]
   
         chooseRow = sample(length(mProb), size=1, prob=mProb)
         tmpMother = mDipTb[chooseRow,]
   
         tmpChild = matrix( c(tmpFather[1], tmpMother[1], tmpFather[1], tmpMother[2],
                              tmpFather[2], tmpMother[1], tmpFather[2], tmpMother[2]), nrow=2, byrow=F)
         t.choice = sample(4, size=1)
   
         trioHapIdx[ss,1:12] =  c(tmpFather, tmpMother, as.vector(tmpChild[, c(t.choice, (1:4)[-t.choice])]))
   
       }
      if(reType){
        return(trioHapIdx)
      }else{
        return(trioHapIdx[,1:12, drop=F])
      }
      
    }else{
      ## no parents is compMiss, child is NOT compMiss, get all and exclude the fDipTb and mDipTb for those doesn't fit cDipTb

      if(ifD) print( "no parents is compMiss, child is NOT compMiss, get all and exclude the fDipTb and mDipTb for those doesn't fit cDipTb")
      if(ifD) print( "Type 6: F(d)M(d)C(d)")
      ## if none is completed missing
      cHapFiltered = procSemiAugMap(appVarNames, cHapBkInfoMap, snpLen)

      cDipTb = fHapBkIdx2DipTb(appVarNames, cHapBkInfoMap, idxList=cHapFiltered, snpCoding,
                                      reqIn=NULL, reqDigits = NULL, expression=cHapBkInfoMap$ori, snpLen)

##       print("###")
##       print(cDipTb)
##       print(fDipTb)
##       print(mDipTb)
      
      ## need to exclude the fDipTb and mDipTb for those doesn't fit cDipTb 
      tmpDip = exParentDip(cDipTb, fDipTb, fHapBkInfoMap,  snpLen)
      cDipTb = tmpDip$childTb
      fDipTb = tmpDip$parTb
#      print(tmpDip)
    
      tmpDip = exParentDip(cDipTb, mDipTb, mHapBkInfoMap, snpLen)
      cDipTb = tmpDip$childTb
      mDipTb = tmpDip$parTb
#      print(tmpDip)
      
      if(!is.null(logF)){
           logl(logF, paste("Possible dip for father with some data: row=", length(fDipTb)/2))
           logl(logF, paste(util.matrix.cat(fDipTb, 1:2, sep="."), collapse="; ", sep=""))
      
           logl(logF, paste("Possible dip for mother with some data: row=", length(mDipTb)/2))
           logl(logF, paste(util.matrix.cat(mDipTb, 1:2, sep="."), collapse="; ", sep=""))
    
           logl(logF, paste("Possible dip for child with some data: row=", length(cDipTb)/2))
           logl(logF, paste(util.matrix.cat(cDipTb, 1:2, sep="."), collapse="; ", sep=""))
      }
    
      ## need to refilter/standandize the pair probability
      ## obtain the diplotype map and probability
      
      prob1 = exDipProbSemiAugMap(fDipTb, semiMapFrame,
                         hapBkOnlyMap.vars$resiProbCol, hapBkOnlyMap.vars$augIdxCol, hapBkOnlyMap.vars$probCol, snpLen)
      
      prob2 = exDipProbSemiAugMap(mDipTb, semiMapFrame,
                         hapBkOnlyMap.vars$resiProbCol, hapBkOnlyMap.vars$augIdxCol, hapBkOnlyMap.vars$probCol, snpLen)
    
      if(ifD) {
        print("Possible diplotype (hap pair) indexes for father (fDipTb); and probability")
        print(fDipTb)
        print(prob1)
        print("Possible diplotype (hap pair) indexes for mother (mDipTb); and probability")
        print(mDipTb)
        print(prob2)
        print("Possible diplotype (hap pair) indexes for child  (cDipTb); ")
        print(cDipTb)
      }

      ## otherwise, build the mating table based on children geno      
      mating6hap = genMatingTBCondOnChild3(appVarNames,  child= cDipTb,
                             par1=fDipTb, par2=mDipTb, prob1, prob2, logF=logF, job=job)

      othChildtt = apply(mating6hap, 1, FUN=find.PsudoControlHap)
      
      trioHapIdx =  cbind(mating6hap[ ,1:4, drop=F],  t(othChildtt))

      #if (ifD) print("@@@@@ return obj@@@@@")
      if(reType){
        #if(ifD) print(  cbind(trioHapIdx, rep(6, job)) )
        return(cbind(trioHapIdx, rep(6, job)))
      }else{
        #if(ifD) print(trioHapIdx)
        return(trioHapIdx)
      }
      
    }
  }

  ## if both parents are completely missing
  if( compMissing.trio[1] & compMissing.trio[2] ){
    if (compMissing.trio[3]){

      if(ifD) print("both parents are compMiss, child is compMiss, sample parents")
      if(ifD) print("Type 1: F(CM)M(CM)C(CM)")

      trioHapIdx = matrix(NA, nrow=job, ncol=13)
      trioHapIdx[,13]=rep(1, job)
      for(ss in 1:job){
        ## both parents are compMiss, child is compMiss, sample parents.
        tmpFather = sampleDipSemiAugMap(semiMapFrame, hapBkOnlyMap.vars$resiProbCol, hapBkOnlyMap.vars$augIdxCol, hapBkOnlyMap.vars$probCol, snpLen)
        tmpMother = sampleDipSemiAugMap(semiMapFrame, hapBkOnlyMap.vars$resiProbCol, hapBkOnlyMap.vars$augIdxCol, hapBkOnlyMap.vars$probCol, snpLen)
  
        tmpChild = matrix( c(tmpFather[1], tmpMother[1], tmpFather[1], tmpMother[2],
                             tmpFather[2], tmpMother[1], tmpFather[2], tmpMother[2]), nrow=2, byrow=F)
        t.choice = sample(1:4, size=1)
  
        trioHapIdx[ss,1:12] =  c(tmpFather, tmpMother, as.vector(tmpChild[, c(t.choice, (1:4)[-t.choice]) ]))
  
      }

      
      if(reType){
        return(trioHapIdx)
      }else{
        return(trioHapIdx[, 1:12, drop=F])
      }
      
    }else{

      if(ifD) print("both parents are compMiss, child is NOT compMiss, sample kids first, then imput parent (ESp method)")
      if(ifD) print("Type 4: F(CM)M(CM)C(d)")
      ## both parents are compMiss, child is NOT compMiss, sample kids first, then imput parent (ESp method)

      supHapProb = getHapProb.semiMapFrame( semiMapFrame, hapBkOnlyMap.vars$resiProbCol,
                           hapBkOnlyMap.vars$augIdxCol, hapBkOnlyMap.vars$probCol, snpLen)

      ## for child not completely missing, need to restandardize the other parents' hap freq
      cHapFiltered = procSemiAugMap(appVarNames, cHapBkInfoMap, snpLen)
      ## need one function to FIX!!!FIX it 
      cDipTb = fHapBkIdx2DipTb(appVarNames, cHapBkInfoMap, idxList=cHapFiltered, snpCoding,
                                  reqIn=NULL, reqDigits = NULL, expression=cHapBkInfoMap$ori, snpLen)

      ## sample child
      cProb = exDipProbSemiAugMap(cDipTb, semiMapFrame=semiMapFrame,
                                  hapBkOnlyMap.vars$resiProbCol, hapBkOnlyMap.vars$augIdxCol,
                                  hapBkOnlyMap.vars$probCol, snpLen)

      trioHapIdx = matrix(NA, nrow=job, ncol=13)
      trioHapIdx[,13]=rep(4, job)
      for(ss in 1:job){
        chooseRow = sample(length(cProb), size=1, prob=cProb)
        hapPair = cDipTb[chooseRow,]
  
        if (hapPair[1]==hapPair[2]){
          matingTbl = ESp.impuParent.homoKid(hapPair[1], supHapProb)
        }else{
          matingTbl = ESp.impuParent.heterKid(hapPair, supHapProb)
        }
  
        trioHapIdxtt = c(matingTbl, hapPair)
        trioHapIdx[ss,1:12] =  c(matingTbl,  find.PsudoControlHap(trioHap6=trioHapIdxtt))
      }

      if(reType){
        return(trioHapIdx)
      }else{
        return(trioHapIdx[,1:12,drop=F])
      }
    }
  }
  

  ## if only one parent is completely missing
  if( sum(compMissing.trio[1:2])==1 ){
    if(ifD) print( "only one parent is completely missing, need to sample the non-comp-missing parents from a restricted list")
    
    ## need to sample the non-comp-missing parents from a restricted list
    nonMparent = trioBlock[!compMissing.trio, ,drop=F][1, ]
    pHapBkInfoMap = hapGenoBlockProc(nonMparent,  snpCoding=snpCoding)
    pHapFiltered = procSemiAugMap(appVarNames, pHapBkInfoMap, snpLen)

#     print(pHapBkInfoMap)
#     print(pHapFiltered)
    ## obtain completed parent filtered dip, knock off the one doesn't match with homozygous index
    pDipTb = fHapBkIdx2DipTb(appVarNames, pHapBkInfoMap,
                                     idxList=pHapFiltered, snpCoding = snpCoding,
                                     reqIn=reqIn, reqDigits = reqDig, expression=pHapBkInfoMap$ori, snpLen)

    if (compMissing.trio[3]){
      ## one parents is compMiss, child is compMiss, sample parents.
      ## just use the popu hap freq for the missing parent
      if(ifD) print( "one parents is compMiss, child is compMiss, sample parents ")
      if(ifD) print("Type 3: F(CM)M(d)C(CM)")
      missingParent = sampleDipSemiAugMap(semiMapFrame, hapBkOnlyMap.vars$resiProbCol, hapBkOnlyMap.vars$augIdxCol, hapBkOnlyMap.vars$probCol, snpLen)

      ## sample the non missing parent
      pProb = exDipProbSemiAugMap(pDipTb, semiMapFrame=semiMapFrame,
                                  hapBkOnlyMap.vars$resiProbCol,
                                  hapBkOnlyMap.vars$augIdxCol,
                                  hapBkOnlyMap.vars$probCol, snpLen)

      trioHapIdx = matrix(NA, nrow=job, ncol=13)
      trioHapIdx[,13]=rep(3, job)
      for(ss in 1:job){
        
        chooseRow = sample(length(pProb), size=1, prob=pProb)
        nomissingParent = pDipTb[chooseRow,]
  
        if(compMissing.trio[1]){
          tmpFather = missingParent
          tmpMother = nomissingParent
        }else{
          tmpFather = nomissingParent
          tmpMother = missingParent
        }
  
        tmpChild = matrix(  c(tmpFather[1], tmpMother[1], tmpFather[1], tmpMother[2],
                              tmpFather[2], tmpMother[1], tmpFather[2], tmpMother[2]), nrow=2, byrow=F)
        
        t.choice = sample(4, size=1)
  
        trioHapIdx[ss,1:12] =  c(tmpFather, tmpMother, as.vector(tmpChild[, c(t.choice, (1:4)[-t.choice]) ]))
      }
      
      if(reType){
        return(trioHapIdx)
      }else{
        return(trioHapIdx[,1:12,drop=F])
      }
      
    }else{

      if(ifD) print( "one parents is compMiss, child is NOT compMiss, use ESp method")
      if(ifD) print("Type 5: F(CM)M(d)C(d)")
      ## one parents is compMiss, child is NOT compMiss, use ESp method.

      ## for child not completely missing, need to restandardize the other parents' hap freq
      cHapFiltered = procSemiAugMap(appVarNames, cHapBkInfoMap, snpLen)
      ## need one function to FIX!!!FIX it 
      cDipTb = fHapBkIdx2DipTb(appVarNames, cHapBkInfoMap, idxList=cHapFiltered, snpCoding,
                                  reqIn=NULL, reqDigits = NULL, expression=cHapBkInfoMap$ori, snpLen)

      #print(cHapFiltered)
      #print(cDipTb) 
      trioHapIdx.f = ESp.impu1Par(
        othParPairs=pDipTb, childPairs=cDipTb,  semiMapFrame,
        resiProbCol=hapBkOnlyMap.vars$resiProbCol,
        augIdxCol=hapBkOnlyMap.vars$augIdxCol,
        probCol=hapBkOnlyMap.vars$probCol,
        snpLen, reType=reType, job=job)

      if(ifD) print(paste("Outcome hap:", paste(trioHapIdx.f, collapse=";")))
      
      ## need to know switch parents if the mom is missing
      if( compMissing.trio[1] ){
        othChildtt = apply(trioHapIdx.f, 1, FUN=find.PsudoControlHap)      
        trioHapIdx =  cbind(trioHapIdx.f[, 1:4, drop=F],  t(othChildtt))
      }else{
        othChildtt = apply(trioHapIdx.f, 1, FUN=find.PsudoControlHap)
        trioHapIdx =  cbind(trioHapIdx.f[, c(3,4, 1,2), drop=F],  t(othChildtt))
      }

      if(reType){
        return(cbind(trioHapIdx, 5+trioHapIdx.f[,7]/10))
      }else{
        return(trioHapIdx)
      }
 
    }
  }
  
}


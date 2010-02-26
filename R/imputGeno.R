imputGeno <-
function(trioBlock, job=1, genoProb, popuProb, data.order, snpCoding){
  fStr = "[imputGeno]:"
  ifD = F
  if(ifD) print(qp(fStr, "start"))
  if( min( c(snpCoding==c(0,1,2,3),data.order=="FMC")) <1 )
     stop (paste(fStr, "Data configuration is not right:\n", "snpCoding=[", paste(snpCoding, collapse=";", sep=""),
                 "]", "trioBlock order=[", data.order, "]", sep=""))

  ##CAUTION, need to flip the popuProb, it is in the sequence for 11, 22, 12
  ## and re-standardize if freq==0
  tmpPop = popuProb/sum(popuProb)
  zeroPop = tmpPop==0
  #print(zeroPop)

  if(sum(zeroPop)>0){
    tmpPop = tmpPop/1.0001
    
    tmpPop[zeroPop]=.0001/sum(zeroPop)

#    if(ifD){
#      print(popuProb)
#      print(tmpPop)
#      print(sum(tmpPop))
#    }
    popuProb=tmpPop
  }
  
  popuProb=popuProb[c(1,3,2)]
  
## first process the block with complete genotypes
  if( sum(trioBlock==snpCoding[1])==0 ){
    seqpar = paste(trioBlock[1], trioBlock[2], sep="-")
    bench = paste(genoProb[,1], genoProb[,2], sep="-")
    matchPos = match(seqpar, bench)

    famRe = matrix(NA, nrow=job, ncol=6)

    
    for( i in 1:job){
      famRe[i, 1:2]= trioBlock[1:2]
      kid = unlist(lapply(1:3, FUN=function(item, repp){
        rep(item, repp[item])
      }, repp=genoProb[matchPos, 3:5]*4))
  
      newlyM = match(trioBlock[3], kid)
  
      allKid = kid[ c(newlyM[1], kid[-newlyM][sample(1:3, size=3, replace=F)])  ]

      famRe[i, 3:6]=allKid      
    }

    return(famRe)
  }

  
  finalIdx = rep(T, nrow(genoProb))
  ## filter out the parents
  if(trioBlock[1]!=snpCoding[1]){
    matched = genoProb[,1]==trioBlock[1]
    finalIdx = finalIdx & matched
  }
  if(ifD) print(finalIdx)
  if(trioBlock[2]!=snpCoding[1]){
    matched = genoProb[,2]==trioBlock[2]
    finalIdx = finalIdx & matched
  }
  if(ifD) print(finalIdx)
  cProb = rep(1, nrow(genoProb))
  if(trioBlock[3]!=snpCoding[1]){
    ## hard code!!! hard code geno coding need to be 1,2,3 to correspond to the child geno
    ## indicated by genoProb
    matched = genoProb[, trioBlock[3]+2]!=0

    cProb = genoProb[, trioBlock[3]+2]
    
    finalIdx = finalIdx & matched
  }
  if(ifD) print(finalIdx)

  if(sum(finalIdx)==0) stop(paste(fStr,
          "Mendelian error! trio geno=[", paste(trioBlock, collapse=".", sep=""),
          "] for order =[", order, "]", sep=""))
#  print("pop")
#  print(popuProb)
  fProb = popuProb[genoProb[,1]]
  mProb = popuProb[genoProb[,2]]
  jProb = fProb*mProb
  jProb = jProb/sum(jProb)

  ## check if estimated genotype return no genotype for this one
  jointProb = sum(jProb[finalIdx])

  if(jointProb==0) {
    jProb[!finalIdx]=0
    jProb = jProb*cProb 
    jProb[finalIdx]=rep(1/sum(finalIdx), sum(finalIdx))
    print("should not happen")
    #print(trioBlock)
    stop("should not happen")
    
  }else{
    jProb[!finalIdx]=0
    jProb = jProb*cProb
    jProb = jProb/sum(jProb)
  }
  
  if(ifD) print(cbind(fProb, mProb, jProb, finalIdx, genoProb[,1:2]))

  allKid = t(apply(genoProb[, 3:5]*4, 1, FUN=function(row){
                  a = rep(c(1, 2, 3), times=row); a}))
  #print(allKid)
  #print(jProb)
  
  re6Geno = matrix(NA, ncol=6, nrow=job)

  for( ss in 1:job){
    
    chooseRow = sample(1:9, size=1, prob = jProb)
  
    re = genoProb[chooseRow, 1:2]
    if(trioBlock[3]!=snpCoding[1]){
      childGeno = trioBlock[3]
    }else{
      ## hard code!!! hard code, geno coding need to be 1,2,3 to correspond to the child geno
      ## indicated by genoProb
      childGeno = sample(1:3, size=1, prob = genoProb[chooseRow, 3:5])
  
    }

    oth = allKid[chooseRow,]
    #print(oth)
    othleft = oth [- which(oth==childGeno)[1] ]
    re6Geno[ss,] = c(re, childGeno, othleft)
  
  }
  #print(re6Geno)
  return(re6Geno)
}


matTBCondOnChild3.finalProc <-
function(maxTbl, counter, maxMateTbl, probVec, prob1, prob2, job, logF){
  #print("Run finalProc")
  ifD = F

  cutoff = min(counter, 20)
  ## in the maxMateTbl, the first column is index within par1 and the second column is index within par2
  if (ifD) print(maxTbl[1:100,])
  if (ifD) print(cbind(maxMateTbl[1:cutoff, ], probVec[1:cutoff]))
  prob1n = rep(0, length(prob1))
  prob1n[unique(maxMateTbl[1:counter,1])]= prob1[unique(maxMateTbl[1:counter,1])]

  prob1n = prob1n/sum(prob1n)
  onePDipProb = prob1n[maxMateTbl[1:counter,1]]
  if(ifD) print(paste("onePDipProb=", paste(round(onePDipProb,3)[1:cutoff], collapse=";", sep="")))


  prob2n = rep(0, length(prob2))
  prob2n[unique(maxMateTbl[1:counter,2])]= prob2[unique(maxMateTbl[1:counter,2])]

  prob2n = prob2n/sum(prob2n)
  othPDipProb = prob2n[maxMateTbl[1:counter,2]]
  if(ifD) print(paste("othPDipProb=", paste(round(othPDipProb,3)[1:cutoff], collapse=";", sep="")))

  prob = onePDipProb * othPDipProb *  probVec[1:counter]
  
  prob = prob/sum(prob)

  if(ifD) print(paste("final prob=", paste(round(prob,3)[1:cutoff], collapse=";", sep="")))

  ## sample the row
  hap6idx = matrix(NA, nrow=job, ncol=6)
  if(ifD) print(hap6idx)
  for(ss in 1:job){
    chooseRow = sample(1:counter, size=1, prob=prob)
    #print(chooseRow)
    hap6idx[ss,] = maxTbl[chooseRow, , drop=F]
    #print(hap6idx)
  }

#   if(!is.null(logF)){
#     if(cutoff==20){
#       bkMapStr = paste(fStr, "\nCutoffed possible dip mating table:\n",
#                      paste(wrComTbl( cbind(maxTbl[1:cutoff, ,drop=F], prob[1:cutoff], probVec[1:cutoff]),
#                               colName=c("f1", "f2", "m1", "m2", "c1", "c2", "prob", "probVec")), collapse="\n", sep=""),
#                      sep="")
#     }else{
#       bkMapStr = paste(fStr, "\nCompleted poss dip mating table:\n",
#                      paste(wrComTbl( cbind(maxTbl[1:cutoff, ,drop=F], prob[1:cutoff], probVec[1:cutoff]),
#                               colName=c("f1", "f2", "m1", "m2", "c1", "c2", "prob", "probVec")), collapse="\n", sep=""),
#                      sep="")
# 
#     }
#     logl(logF, bkMapStr)
#     logl(logF, paste(fStr, "choose row num=", chooseRow, sep=""))
#     #logl(logF, paste("Parents index=[", paste(as.vector(hap6idx), collapse=".", sep=""), "]", sep=""))
#   }

  if(ifD) {
    tttt = cbind(maxTbl[1:cutoff, ,drop=F], prob[1:cutoff])
    print(tttt)
    print(hap6idx)
  }
  rm(maxTbl)
  rm(maxMateTbl)
  rm(probVec) 

  return( hap6idx )

}


bkMap.HRCB.famMap <-
function(bkMapS, rule, newColName,  ifS="simuDirectInfo", baseName=NULL){
     ifD = F

     ct.shap = cumprod(bkMapS$bkLens)[length(bkMapS$bkLens)]
     ct.sdip = .5*ct.shap*(1+ct.shap)
     ct.mrow = .5*ct.sdip*(1+ct.sdip)

     if(ct.mrow > 5*10^6) warning(qp("Total number of super-haplotype exceed 5*10^6. It is likely to exceed the memory limit."))

     if(ct.mrow > 5*10^7) stop(qp("Total number of super-haplotype exceed 5*10^7. It exceeds the memory limit."))

     linkedHap = bkMap.superHRCB(bkMapS)
     hapExp = apply(linkedHap$linkBkExp, 1, paste, collapse="")
     prob = linkedHap$prob/sum(linkedHap$prob)
   
     hapCt = length(prob)
     
     ab = hap2GenoBlock(hapExp, prob, snpCt = bkMapS$snpLen, alleleCode = bkMapS$alleleCode)

     risk = HRCB.applyRule(ab$tb$genoExp, rule, newColName)

     ## ruleMet only return the applied result for the last rule in the list
     ## associate the hap pair index with prob is used for the kids
     if(length(rule$slist)==1){
       tb = cbind(ab$tb, ruleMet=risk$treePred, risk.Prob=risk$riskProb)
     }else{
       tb = data.frame(ab$tb, risk.Prob=risk$riskProb)
     }

     if(ifD) print(str(tb))
     if(ifD) print(tb)

     #print("test")

     ## generat kids hap ids pair
     rowCt = nrow(tb)

     if(ifD) print(rowCt)
     matRowCt = (rowCt^2-rowCt)/2+rowCt
     matingRowIdxCorn = util.it.smallLargeIdx(rowCt, keep.same=F)
     if(ifD) print(str(matingRowIdxCorn))
     matingRowIdxDiag = matrix(rep(1:rowCt, each=2), ncol=2, byrow=T)
     if(ifD) print(str(matingRowIdxDiag))

     matingRowIdx = rbind(matingRowIdxCorn, matingRowIdxDiag)
     if(ifD) print(matingRowIdx[1:20,])

     tmpProb = tb$prob
     matingPCor = matrix(tmpProb[matingRowIdxCorn], ncol=2, byrow=F)
     matingPCorn = 2*matingPCor[,1]*matingPCor[,2]
     matingPDiag = tmpProb^2
     matingP = c(matingPCorn, matingPDiag)
     if(ifD) print(matingP[1:20])
     
     ## HARD CODE!!! NEED to be blocked
##****  subs = 1:10
##****  matingRowIdx = matingRowIdx[subs,]
##****  matingP = matingP[subs]

     tmpTb = matrix(unlist(tb[,c(1,2)]), ncol=2, byrow=F)
     if(ifD) print(str(tmpTb))
    
     fa.hapIdx = tmpTb[matingRowIdx[,1], ]
     ma.hapIdx = tmpTb[matingRowIdx[,2], ]

     if(ifD) {
       print("Sample parents hap idxes:")
       print(fa.hapIdx[1:10,])
       print(ma.hapIdx[1:10,])
     }

     fam.map = matrix(NA, ncol=12, nrow=matRowCt )
     fam.map[,c(1,2)]=fa.hapIdx
     fam.map[,c(3,4)]=ma.hapIdx
     fam.map[,5]=pmin( fa.hapIdx[,1], ma.hapIdx[,1])
     fam.map[,6]=pmax( fa.hapIdx[,1], ma.hapIdx[,1])
     fam.map[,7]=pmin( fa.hapIdx[,1], ma.hapIdx[,2])
     fam.map[,8]=pmax( fa.hapIdx[,1], ma.hapIdx[,2])
     fam.map[,9]=pmin( fa.hapIdx[,2], ma.hapIdx[,1])
     fam.map[,10]=pmax( fa.hapIdx[,2], ma.hapIdx[,1])
     fam.map[,11]=pmin( fa.hapIdx[,2], ma.hapIdx[,2])
     fam.map[,12]=pmax( fa.hapIdx[,2], ma.hapIdx[,2])

     if(ifD) print(fam.map[1:10,])

     ## obtain the disease prob given the hap idx
     kids.hapIdx = matrix(t(fam.map[,5:12]), nrow=2, byrow=F, ncol=4*matRowCt)
     if(ifD) print(kids.hapIdx[,1:10])
     kids.hapIdx = t(kids.hapIdx)
     if(ifD) print(kids.hapIdx[1:10,])
     if(ifD) print(str(kids.hapIdx))

### instead of using the inefficient matching, we will use the virtual mapping function 
#     kids.matchedRow = unlist(apply(kids.hapIdx[1:20,], 1, util.vec.matchVecIdx, vec=t( tmpTb ),  vecLen=rowCt*2, #benchLen=2))

     kids.matchedRow = unlist(apply(kids.hapIdx, 1, util.it.triMatch, len=hapCt))

     if(ifD) print(kids.matchedRow[1:20])
     # kids.p = matrix(  rep(.25*matingP, each=4), ncol=4, byrow=T)
     # print( paste("Check:: kids p sum (before standardization)=", sum(kids.p)))
     kids.p = matingP/sum(matingP)
  
     ## check
     #if(ifD)  print( paste("Check:: kids p sum=", sum(kids.p)))
     #fam.map = data.frame(fa.hapIdx, ma.hapIdx, ch1.h, ch2.h, ch3.h, ch4.h, rowProb = matingP, kids.p)
     #if(ifD) print(fam.map[1:10,])
  
     kids.risk = matrix( risk$riskProb [kids.matchedRow], ncol=4, byrow=T)
     kids.matchedRow = matrix(kids.matchedRow, ncol=4, byrow=T)
     if(ifD) {
       print("kids risk")
       print(str(kids.risk))
     }
     kids.risk[,1] = (.25*kids.p)* kids.risk[,1]
     kids.risk[,2] = (.25*kids.p)* kids.risk[,2]
     kids.risk[,3] = (.25*kids.p)* kids.risk[,3]
     kids.risk[,4] = (.25*kids.p)* kids.risk[,4]
  
     ## check
     #print( paste("Check:: kids risk sum=", sum(kids.risk)))
     
     fam.map = cbind(fam.map, matingP, kids.risk)
     colnames(fam.map) = c("f.hap1", "f.hap2", "m.hap1", "m.hap2",
                         "c1.hap1", "c1.hap2", "c2.hap1", "c2.hap2",
                         "c3.hap1", "c3.hap2", "c4.hap1", "c4.hap2",
                         "rowProb", 
                         "c1Risk", "c2Risk", "c3Risk", "c4Risk")
  
     if(ifD) print(round(fam.map[1:10,],5))

     if(!is.null(ifS)) {
       
       #write.csv(tb, file=paste(ifS, "hap2geno.csv", sep=""))
       #write.csv(fam.map,  file=paste(ifS, "hapMating.csv", sep=""))
     }
  
     matingTbInfo = list(matRowCt=matRowCt, kids.risk=kids.risk,
                 matingRowIdx = matingRowIdx,
                 hap2genoMap = ab$tb,
                 snpLen = bkMapS$snpLen,
                 kids.matchedRow = kids.matchedRow)
     #print(str(matingTbInfo))

     if(!is.null(baseName))  save( matingTbInfo, file=paste(baseName, ".RData", sep=""))

     return( matingTbInfo )
}


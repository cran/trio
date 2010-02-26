HRCB.famMap.spTrio <-
function(caseNo, matingTbName, ifS="simuDirectInfo", reControl =F ){
     ifD = F

     FN = "HRCB.famMap.spTrio" 
     #print(FN)
     matingTbInfo  = get(matingTbName)
     #print(str(matingTbInfo))

     
     matRowCt = matingTbInfo$matRowCt
     #print(qp("matRowCt=", matRowCt*4))
     #print(str(matingTbInfo))
     
     kids.risk = matingTbInfo$kids.risk
     matingRowIdx = matingTbInfo$matingRowIdx
     genoMap = matingTbInfo$hap2genoMap
     snpLen = matingTbInfo$snpLen
     kids.matchedRow = matingTbInfo$kids.matchedRow

     ## sample row
     ##caseNo = 5
     kid.idx.sample = sample( 1:(matRowCt*4), size = caseNo, prob=kids.risk, replace=T)

     ## find the matrix idx for the kids
     row.sample = kid.idx.sample%%matRowCt
     row.sample[row.sample==0]=matRowCt
     case.idx.matchedRow = kids.matchedRow[kid.idx.sample]
     
     if(reControl){
       ## need to calculate the control id
       ## recreat the child id for sampled rows
       controlIdx = unlist(lapply(1:caseNo, FUN=function(i, totalRow, rowIdx, cIdx){
             childIDX = rep(rowIdx[i], 4)+ (0:3)*totalRow
             leftControl = childIDX[ childIDX!=cIdx[i] ]
           }, totalRow = matRowCt, rowIdx=row.sample, cIdx=kid.idx.sample))

       tt.newIdxSeq = cbind( case.idx.matchedRow, matrix( kids.matchedRow[controlIdx], ncol=3, byrow=T))
       genoOth = genoMap[,6][ tt.newIdxSeq  ]
       geno.CC = genoOth[t( matrix(1:(caseNo*4), ncol=4, byrow=F) )   ]
       geno.CCMa = t(sapply(geno.CC, FUN=geno.2dStr2BinaMa, subjectCt=caseNo*4, snpLen=snpLen))
       
     }


     ## generate parents, kids, pesudo controls
     fa.idx.matchedRow = matingRowIdx[,1][row.sample]
     ma.idx.matchedRow = matingRowIdx[,2][row.sample]
     

     if(!is.null(ifS)) {
       sample.idx = data.frame(kid.idx.sample, row.sample, fa.idx.matchedRow, ma.idx.matchedRow, case.idx.matchedRow)
       colnames(sample.idx) = c("idx_kids", "rowIdx_matingTb", "fRowIdx_genoTb",
                                 "mRowIdx_genoTb", "cRowIdx_genoTb")
       write.csv(sample.idx,  file=paste(ifS, "supHap.csv", sep=""))
     }

     ## convert string into digits
     hap.str1 = genoMap[,3][ c(fa.idx.matchedRow, ma.idx.matchedRow, case.idx.matchedRow) ]
     hap.str2 = genoMap[,4][ c(fa.idx.matchedRow, ma.idx.matchedRow, case.idx.matchedRow) ]

     geno = covDipStr2CodedGeno(hap.str1, hap.str2, subjectCt=caseNo*3, snpLen=snpLen, snpCoding=0:3, snpBase=c(0, 1, 2))

     #geno = genoMap[,6][c(fa.idx.matchedRow, ma.idx.matchedRow, case.idx.matchedRow)]

     geno.FMCMa = geno[ t( matrix(1:(caseNo*3), ncol=3, byrow=F) ),   ]
     if(ifD) print(matrix(geno.FMCMa, ncol=1)[1:10,])
     ## convert string into digits
     # geno.FMCMa = t(sapply(geno.FMC, FUN=geno.2dStr2BinaMa, subjectCt=caseNo*3, snpLen=snpLen))
   
     if(!reControl){
        return(geno.FMCMa)
     }

     
    stop ("NOT implemented yet")
    return(list(cc=geno.CCMa, geno.FMCMa=geno.FMCMa))

}


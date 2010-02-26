bkMap.LRCB.spTrio <-
function(bkMap, caseNo, ifS = NULL, reControl=F){

   ifD = F
   # first sample parents
   f.samples = bkMap.spHap(bkMap, caseNo)
   f.samples2 = bkMap.spHap(bkMap, caseNo)

   m.samples = bkMap.spHap(bkMap, caseNo)
   m.samples2 = bkMap.spHap(bkMap, caseNo)

   #dim(par)
   fa = covDipStr2CodedGeno(f.samples$subjects, f.samples2$subjects, subjectCt=caseNo, snpLen=bkMap$snpLen, snpCoding=0:3, snpBase=c(0, bkMap$alleleCode))
   ma = covDipStr2CodedGeno(m.samples$subjects, m.samples2$subjects, subjectCt=caseNo, snpLen=bkMap$snpLen, snpCoding=0:3, snpBase=c(0, bkMap$alleleCode))

   ## shuffle the 
   par = cbind(as.vector(f.samples$subjects), as.vector(f.samples2$subjects),
               as.vector(m.samples$subjects), as.vector(m.samples2$subjects))

   parIn =  cbind(as.vector(f.samples$subjectsIn), as.vector(f.samples2$subjectsIn),
               as.vector(m.samples$subjectsIn), as.vector(m.samples2$subjectsIn))

   #print(paste("!reControl: par dim:", paste(dim(par), collapse=" by ", sep="")))

   
   ## shuffle the random kid
   ## get index, !!!take the same happair combination for all blocks and all trio
   baseIdx =  matrix(c(1,3,2,3,1,4,2,4), nrow=2, byrow=F)

   randomRowIdx = sample(1:4, size=1, replace=F)
   newchild = baseIdx[,randomRowIdx]

   chExp = par[, newchild]
   chIn = parIn[, newchild]

   #print(paste("!reControl: chExp dim:", paste(dim(chExp), collapse=" by ", sep="")))

   #print(dim(chExp))
   #print(dim(chIn))
   #print(dim(newchild))

   child1 = matrix(chExp[,1], nrow=caseNo, byrow=F)
   child2 = matrix(chExp[,2], nrow=caseNo, byrow=F)

   affChild = covDipStr2CodedGeno(child1, child2, subjectCt=caseNo, snpLen=bkMap$snpLen, snpCoding=0:3, snpBase=c(0, bkMap$alleleCode))
   trioSetts = rbind(fa, ma, affChild)
   #print(paste("!reControl: trioSetts dim:", paste(dim(trioSetts), collapse=" by ", sep="")))

   if(!reControl){     
     rowNewSeq = t(matrix(1:(caseNo*3), ncol=3, nrow=caseNo, byrow=F))
     geno.FMCMa = trioSetts[rowNewSeq,]
   }else{

     othRowIdx = (1:4)[ -randomRowIdx ]
     
     exHapIdx = baseIdx[,othRowIdx]

     tt = F
     if(tt) print(paste("othRowIdx=", paste(othRowIdx, collapse=";")))
     if(tt) print(paste("expHapIdx=", paste(exHapIdx, collapse=";")))

     for( i in 1:3){
          childExpIdx = exHapIdx[,i]
          #print(childExpIdx)
          othChildExp = par[,childExpIdx]
          #print(paste("!reControl: othChildExp dim:", paste(dim(othChildExp), collapse=" by ", sep="")))
          childOth1 = matrix(othChildExp[,1], nrow=caseNo, byrow=F)
          childOth2 = matrix(othChildExp[,2], nrow=caseNo, byrow=F)

          affChild =
            covDipStr2CodedGeno(childOth1, childOth2, subjectCt=caseNo, snpLen=bkMap$snpLen, snpCoding=0:3, snpBase=c(0, bkMap$alleleCode))
          trioSetts = rbind(trioSetts, affChild)
     }
     #print(paste("!reControl: trioSetts dim:", paste(dim(trioSetts), collapse=" by ", sep="")))
     
     rowNewSeq = t(matrix(1:(caseNo*6), ncol=6, nrow=caseNo, byrow=F))
     #print(trioSetts)
     geno.FMCMa = trioSetts[rowNewSeq,]
     

   }

     ## rearrange the data, so the three subject in a family goes together.
   if(!is.null(ifS)){
       ## retain trio Index
       childIn1 = matrix(chIn[,1], nrow=caseNo, byrow=F)
       childIn2 = matrix(chIn[,2], nrow=caseNo, byrow=F)
       faIn1 = matrix(parIn[,1],  nrow=caseNo, byrow=F)
       faIn2 = matrix(parIn[,2],  nrow=caseNo, byrow=F)
       maIn1 = matrix(parIn[,3],  nrow=caseNo, byrow=F)
       maIn2 = matrix(parIn[,4],  nrow=caseNo, byrow=F)
  
       childIn = util.matrix.col.shuffle2(childIn1, childIn2)
       faIn = util.matrix.col.shuffle2(faIn1, faIn2)
       maIn = util.matrix.col.shuffle2(maIn1, maIn2)
  
       allIn = rbind(faIn, maIn, childIn)

       rowNewSeq = t(matrix(1:(caseNo*3), ncol=3, nrow=caseNo, byrow=F))
       allIn = allIn[rowNewSeq,]

       write.table(allIn,  file=paste(ifS, "supHap.csv", sep=""), col.names=F, row.names=F, sep=",")         
   }
   return(geno.FMCMa)
   
#      ## not implemented
#      ## shuffle the random kid
#      child = apply(par, 1, FUN=function(row){
#        idx = matrix(c(1,3,2,3,1,4,2,4), nrow=2, byrow=F)
#        randomChildIdx = sample(1:4, size=4, replace=F)
#        idx = idx[,randomChildIdx]
#   
#        newchild = row[idx]
#      })
#      child1.1 = matrix(child[1,], nrow=caseNo, byrow=F)
#      child1.2 = matrix(child[2,], nrow=caseNo, byrow=F)
# 
#      child2.1 = matrix(child[3,], nrow=caseNo, byrow=F)
#      child2.2 = matrix(child[4,], nrow=caseNo, byrow=F)
# 
#      child3.1 = matrix(child[5,], nrow=caseNo, byrow=F)
#      child3.2 = matrix(child[6,], nrow=caseNo, byrow=F)
# 
#      child4.1 = matrix(child[7,], nrow=caseNo, byrow=F)
#      child4.2 = matrix(child[8,], nrow=caseNo, byrow=F)
# 
#      affChild1 = covDipStr2CodedGeno(child1.1, child1.2, subjectCt=caseNo, snpLen=bkMap$snpLen, snpCoding=0:3, snpBase=c(0, bkMap$alleleCode))
#      affChild2 = covDipStr2CodedGeno(child2.1, child2.2, subjectCt=caseNo, snpLen=bkMap$snpLen, snpCoding=0:3, snpBase=c(0, bkMap$alleleCode))
#      affChild3 = covDipStr2CodedGeno(child3.1, child3.2, subjectCt=caseNo, snpLen=bkMap$snpLen, snpCoding=0:3, snpBase=c(0, bkMap$alleleCode))
#      affChild4 = covDipStr2CodedGeno(child4.1, child4.2, subjectCt=caseNo, snpLen=bkMap$snpLen, snpCoding=0:3, snpBase=c(0, bkMap$alleleCode))
#   
#      allChild = rbind(affChild1, affChild2, affChild3, affChild4)
# 
#      # child is SNPct*4 by caseNo matrix
#      newSeq = matrix(1:(4*caseNo), ncol=caseNo, byrow=T)
#      newChild = allChild[newSeq,]
# 
#      trioSetts = rbind(fa, ma, affChild1)
#      rowNewSeq = t(matrix(1:(caseNo*3), ncol=3, nrow=caseNo, byrow=F))
# 
#      geno.FMCMa = trioSetts[rowNewSeq,]
#      return(trio=list(geno.FMCMa=geno.FMCMa, cc = newChild))

}


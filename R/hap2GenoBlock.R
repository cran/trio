hap2GenoBlock <-
function(hapExp = c("11", "12", "21", "22"), hapProb = rep(1/length(hapExp), length(hapExp)), snpCt = nchar(hapExp[1]), alleleCode = c(1, 2)){

     ## first standardize??
     
     ct = length(hapExp)
#      hapIdxCorner=NULL
#      ## 2 hap combination to form the diplotype, identify each variation by its index
#      for ( i in 1:ct ){
#        for ( j in 1:ct){
#          if(i<j) hapIdxCorner = rbind(hapIdxCorner, c(hap1=i, hap2=j))
#        }
#      }
#      
# 
     ifD = F
     if(ifD) print(paste("hapLen:", ct))
     hapIdxCorner = util.it.smallLargeIdx(ct, keep.same=F)
     hapIdxDiag = matrix(rep(1:ct, times=2), ncol=2, byrow=F)
  
     probCorner = matrix(as.vector(hapProb)[as.vector(hapIdxCorner)], ncol=2, byrow=F)
     probDiag = matrix(as.vector(hapProb)[as.vector(hapIdxDiag)], ncol=2, byrow=F)
     joinCorner = 2*probCorner[,1]*probCorner[,2]
     joinDiag = probDiag[,1]*probDiag[,2]
     
     expCorner = matrix(hapExp[hapIdxCorner], ncol=2, byrow=F)
     expDiag = matrix(hapExp[hapIdxDiag], ncol=2, byrow=F)
     
     hapIdx = rbind(hapIdxCorner, hapIdxDiag)
     colnames(hapIdx) = c("id1", "id2")
     hapExp = rbind(expCorner, expDiag)
     colnames(hapExp) = c("hExp1", "hExp2")

     if(ifD) print(hapExp)
     # get the one-digit coding
     # when add to digit, get c(2,3,4)-1=c(1,2,3) the common SNP coding
     # HARD CODE!!!HARD CODE
     bina1 = hapBk2AlleleSeq(hapExp[,1], subjectCt=dim(hapExp)[1], snpLen=snpCt, markdownOne = F)
     bina2 = hapBk2AlleleSeq(hapExp[,2], subjectCt=dim(hapExp)[1], snpLen=snpCt, markdownOne = F)
      if(ifD) print(bina1)
     
     ##  20Dec08Change!!!: straighten coding: output must be of digit 1, 2, 3 for 1-digitCoding, and 1, 2 for 2-digit
     a1 = sum(alleleCode * c(2,0)) # 11-> to 1
     a2 = sum(alleleCode * c(1,1)) # 12-> to 3
     a3 = sum(alleleCode * c(0,2)) # 22-> to 2

     snp1d = bina1+bina2
     snp1d.f = factor(snp1d, levels=c(a1, a3, a2))

     if(ifD) {
       print(c(a1, a2, a3))
       print(";;")
       print(snp1d)
       print(snp1d.f)
       print(is.na(snp1d.f))
     }
     
     if( max(is.na(snp1d.f))==1) stop("Input allelCode does not match with other input.")
     
     snp1d = as.integer(snp1d.f)
     snp1d = matrix(snp1d, ncol=ncol(bina1), nrow=nrow(bina1))

     # get the two-digit coding
     bina1.f = factor(bina1, levels=alleleCode)
     if( max(is.na(bina1.f))==1) stop("Input allelCode does not match with other input.")
     bina1.f = as.integer(bina1.f)

     bina2.f = factor(bina2, levels=alleleCode)
     if( max(is.na(bina2.f))==1) stop("Input allelCode does not match with other input.")
     bina2.f = as.integer(bina2.f)

     bina1 = matrix(bina1.f, ncol=ncol(bina1), nrow=nrow(bina1))
     bina2 = matrix(bina2.f, ncol=ncol(bina1), nrow=nrow(bina1))
     
     binaNew1 = pmin(bina1, bina2)
     binaNew2 = pmax(bina1, bina2)

     if(ifD) print(binaNew1)
     
     genoTypes.m = util.matrix.col.shuffle2(binaNew1, binaNew2)
     genoExp = util.matrix.cat(genoTypes.m, 1:(snpCt*2), sep="")
      
     prob = matrix(c(joinCorner, joinDiag), ncol=1)

     tb = data.frame(hapIdx, hapExp, prob, genoExp)
     
     return(list(tb=tb, genoMa=genoTypes.m, geno1d = snp1d))
}


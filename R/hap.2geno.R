hap.2geno <-
function(hapExp, hapProb, snpBIdx=NA){

    if(is.na(snpBIdx)) snpBIdx = 0
    bina1 = hapBk2AlleleSeq(subjects=hapExp,
                            subjectCt=length(hapExp),
                            snpLen=nchar(hapExp[1]), markdownOne = F)
    
    geno1dFreq = lapply(1:ncol(bina1), FUN=function(col, prob, snpMa){
      snpF = factor(snpMa[,col],levels=1:2)
      
      gFreq = tapply( prob, snpF, sum)
      #print(gFreq)
      gFreq[ is.na(gFreq) ] = 0
      gFreq1 = c(gFreq, 2*gFreq[1])
      gFreq2 = c(gFreq, gFreq[2])
      re = gFreq1*gFreq2
      #print(re)
    }, prob=hapProb  , snpMa = bina1)
    
    genoFreq = util.list.2matrix(list=geno1dFreq, byRow = T)
    colnames(genoFreq)=c("11","22","12")
    rownames(genoFreq)=paste("snp", snpBIdx+1:ncol(bina1), sep="")
    return(genoFreq)
}


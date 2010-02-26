getBackParentGeno <-
function(trioDf, famCol=1, memCol=2, snpIdx=NULL, prefix=NULL, re.child=F){

    if(is.null(snpIdx)){
      snpIdx = (1:ncol(trioDf))[-c(famCol, memCol)]
    }
      ## HARD CODE!!!HARD CODE: know the last three digits are for covariate values
     rowCt = nrow(trioDf)
     trioSNP=trioDf[, c(famCol, memCol, snpIdx)]

     if(re.child) {
       selMa = seq.int(from=3, to=rowCt, by = 3)

     }else{
       ## we know the first two rows in every three rows are parents
       selSeq = seq.int(from=1, to=rowCt, by = 3)
       selSeq2 = selSeq + 1
       selMa = matrix(c(selSeq, selSeq2), ncol=2, byrow = F)
       selMa = as.vector(t(selMa))
     }

     #print(length(selMa))
     trioSNP.id = paste(trioSNP[selMa, 1], trioSNP[selMa, 2], sep="-")
     unique.par = unique(trioSNP.id)

     unique.par.idx = match(unique.par, trioSNP.id)

     #print(length(unique.par))
     
     if(is.null(prefix)){
       if(re.child) {
         re = trioSNP[unique.par.idx,]
       }else{
         re = trioSNP[unique.par.idx,]
       }
       return(re)
     }
     
     if(re.child) {
       write.table(trioSNP[unique.par.idx,], file=paste(prefix, "_child.txt", sep=""), sep=" ",
                   append = FALSE, row.names = F, col.names = F)
     }else{
       write.table(trioSNP[unique.par.idx,], file=paste(prefix, "_parent.txt", sep=""), sep=" ",
                   append = FALSE, row.names = F, col.names = F)
     }
     return(NULL)
}


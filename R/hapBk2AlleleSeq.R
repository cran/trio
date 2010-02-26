hapBk2AlleleSeq <-
function(subjects, subjectCt, snpLen, markdownOne=T){

  if(is.null(dim(subjects))){
     subjects.snp = lapply(subjects, FUN=util.str.2CharArray,  len=snpLen)    
  }else{     
     bkCt = dim(subjects)[2]
     
     subjects.hap = util.matrix.cat(subjects, 1:bkCt, sep="")
     subjects.snp = lapply(subjects.hap, FUN=util.str.2CharArray,  len=snpLen)
   }
  
  if(markdownOne){
    subjects.snp = matrix(as.numeric(unlist(subjects.snp))-1, ncol = snpLen, byrow = T)
  }else{
    subjects.snp = matrix(as.numeric(unlist(subjects.snp)), ncol = snpLen, byrow = T)
  }
  return(subjects.snp)
  
}


geno.2dStr2BinaMa <-
function(expStr, subjectCt=length(expStr), snpLen){

  samplesBina = hapBk2AlleleSeq(subjects=expStr, subjectCt=subjectCt, snpLen=snpLen*2)
  bina1= samplesBina[, seq.int(from=1,to=snpLen*2, by=2), drop=F]
  bina2 = samplesBina[, seq.int(from=2,to=snpLen*2, by=2), drop=F]
   binaNew1 = pmin(bina1, bina2)
   binaNew2 = pmax(bina1, bina2)
   bina = util.matrix.col.shuffle2(binaNew1, binaNew2)

   if(subjectCt==1) bina=matrix(bina, nrow=1)
   return(bina)  

}


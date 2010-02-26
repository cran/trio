bkMap.spGenoSeq <-
function(bkMap, subjectCt){

   # first set of hap
   samples = bkMap.spHap(bkMap, subjectCt)
   samples2 = bkMap.spHap(bkMap, subjectCt)
   if(F) print(cbind(samples[[1]], samples2[[1]]))

   # change hap to genotype
   samplesBina = hapBk2AlleleSeq(subjects=samples$subjects, subjectCt, snpLen=bkMap$snpLen)
   samplesBina2 = hapBk2AlleleSeq(subjects=samples2$subjects, subjectCt, snpLen=bkMap$snpLen)
  
   ## construct the genotypic data
   binaNew1 = pmin(samplesBina, samplesBina2)
   binaNew2 = pmax(samplesBina, samplesBina2)
   bina = util.matrix.col.shuffle2(binaNew1, binaNew2)

   return(bina)
}


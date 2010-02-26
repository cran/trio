hap2genotype.m <-
function(hapMa1, hapMa2, subjectCt, snpLen){
  bina1 = hapBk2AlleleSeq(hapMa1, subjectCt, snpLen, markdownOne = F)
  bina2 = hapBk2AlleleSeq(hapMa2, subjectCt, snpLen, markdownOne = F)
  
  binaNew1 = pmin(bina1, bina2)
  binaNew2 = pmax(bina1, bina2)
  re = util.matrix.col.shuffle2(binaNew1, binaNew2)

  return(re)
}


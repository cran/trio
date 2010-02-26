covDipStr2CodedGeno <-
function(dipStr1, dipStr2, subjectCt, snpLen, snpCoding=c(0,1,2,3), snpBase=c(0,1,2)){
  dip1 =  hapBk2AlleleSeq(dipStr1, subjectCt, snpLen, markdownOne = F)
  dip2 =  hapBk2AlleleSeq(dipStr2, subjectCt, snpLen, markdownOne = F)
  re = covDipBinaMa2CodedGeno(dip1, dip2, subjectCt, snpCoding=snpCoding, snpBase=snpBase)
  return(re)
}


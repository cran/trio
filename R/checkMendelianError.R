checkMendelianError <-
function(codedSNPTrio, snpCoding=c(0,1,2,3)){

  # is the child homo with 1?
  if(codedSNPTrio[3]==snpCoding[2]){
    onePHaveNone = F
    if(codedSNPTrio[1]==snpCoding[3]) onePHaveNone = T
    othPHaveNone = F
    if(codedSNPTrio[2]==snpCoding[3]) othPHaveNone = T
    if(onePHaveNone | othPHaveNone)
      stop("Medelian error for homozygous (1) child with at least one parent homozygous (2)")
  # is the child homo with 2?
  }else if(codedSNPTrio[3]==snpCoding[3]){
    onePHaveNone = F
    if(codedSNPTrio[1]==snpCoding[2]) onePHaveNone = T
    othPHaveNone = F
    if(codedSNPTrio[2]==snpCoding[2]) othPHaveNone = T
    if(onePHaveNone | othPHaveNone)
      stop("Medelian error for homozygous (2) child with at least one parent homozygous (1)")
  # is the child hetero
  }else if(codedSNPTrio[3]==snpCoding[4]){

    if(codedSNPTrio[1]==snpCoding[2]){
      if(codedSNPTrio[2]==snpCoding[2]){
        stop("Medelian error for heterozygous child with two parents homozygous (1)") 
      }
    }else if(codedSNPTrio[1]==snpCoding[3]){
      if(codedSNPTrio[2]==snpCoding[3]){
        stop("Medelian error for heterozygous child with two parents homozygous (2)") 
      }
    }
  # is the child missing
  }

  return(NULL)
}


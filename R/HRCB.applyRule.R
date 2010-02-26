HRCB.applyRule <-
function(subjectsExpMa, rule, col.shuffled){
  ifD = F
  len = length(rule)
  if(is.null(dim(subjectsExpMa))){
    nsub = length(subjectsExpMa)
  }else{
    nsub = nrow(subjectsExpMa)
  }
  varCt = length(col.shuffled)

  if(ifD) print("subjectsExpMa:")
  if(ifD) print(subjectsExpMa[1:10])
  subjectListMa = hapBk2AlleleSeq(subjectsExpMa, nsub, varCt, markdownOne=T)
  print(str(subjectListMa))

  subjectList = util.matrix.2list(subjectListMa)
  
  linear = unlist(rep(rule$inter, nsub))
  slist = rule$slist
  
  for ( i in 1:length(slist)){
    curRule = slist[[i]]
    treePred = binaTree.apply(binaTree=curRule$signal, data=subjectListMa, colnames=col.shuffled, re.final=T)
    #print(treePred)
    linear = linear + treePred * unlist(rep(curRule$coef, nsub))
  }

  if(ifD){
    aaaa = cbind(subjectsExpMa, subjectListMa, treePred )
    write.csv(aaaa, file="HRCB.applyRule.diagnostic.csv")

  }
  #print(linear)
  riskProb =  getP.logodds(linear)
  attributes(riskProb)=NULL
  #print(riskProb)
  return(list(riskProb=riskProb, treePred=treePred))
}


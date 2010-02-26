hap2genotype.bk <-
function(bkDframe, chCol, blockCol, keyCol, expCol, probCol, hapLenCol){
  hapExp = bkDframe[,expCol]
  prob = bkDframe[,probCol]
  dipExp =expand.grid(hapExp, hapExp)
  probExp = expand.grid(prob, prob)
  varCt =  bkDframe[1, hapLenCol]
  genotype = hap2genotype.m(as.character(dipExp[,1]), as.character(dipExp[,2]), nrow(dipExp), varCt)
  probExp = probExp[,1]*probExp[,2]
  geno = util.matrix.cat( genotype, 1:(2*varCt))
  collaps = tapply(probExp, geno, sum)

  nrow = length(collaps)
  
  re = data.frame(ch=rep(bkDframe[1,chCol], nrow), block=rep(bkDframe[1,blockCol], nrow),
                  exp=names(collaps), varCt=rep(2*varCt, nrow), prob=collaps,
                  key=rep(bkDframe[1,keyCol], nrow))
 ## re = list(probList=collaps, varCt = varCt)
  return(re)
}


find.PsudoControlHap <-
function(trioHap6){
  
  child.sort = sort(trioHap6[5:6])
  
  othChild = matrix(trioHap6[1:4][c(1, 3, 1, 4, 2, 3, 2, 4)], ncol=2, byrow=T)
  othChild = apply(othChild, 1, sort)
  
  child.m = NULL
  for( i in 1:4){
    oth = othChild[,i]
    if(sum(child.sort==oth)==2){
      child.m = c(child.m, i)
    }
  }
  if(length(child.m)<1) stop("No match for the child")

  t.choice = child.m[sample(length(child.m), size=1)]
 
  re = othChild[, c(t.choice, (1:4)[-t.choice]) ]

  return(re)
}


genGenoProb <-
function(){
   genoProb = matrix(c(1,1, 1, 0, 0,
                       2,1, 0, 0, 1,
                       3,1,.5, 0,.5,
                       1,2, 0, 0, 1,
                       2,2, 0, 1, 0,
                       3,2, 0,.5,.5,
                       1,3,.5, 0,.5,
                       2,3, 0,.5,.5,
                       3,3,.25,.25, .5), ncol=5, byrow=T)
   ## test
   ## apply(genoProb[,3:5], 1, sum)
   ## apply(genoProb[,3:5], 2, sum)
   return (genoProb)
}


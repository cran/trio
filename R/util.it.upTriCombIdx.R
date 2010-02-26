util.it.upTriCombIdx <-
function(idx, diag=T, re.ordered = F){
     ct = length(idx)
     reRow = (ct^2 - ct)/2

     if(diag){
       reRow = reRow + ct
       re = matrix(NA, ncol=2, nrow=reRow)
       if(re.ordered){
          it = 0
          for ( i in 1:ct ){
            for ( j in i:ct){
                it = it + 1
                re[it,] = range(c(idx[i], idx[j]))
            }
          }
       }else{
          it = 0
          for ( i in 1:ct ){
            for ( j in i:ct){
                it = it + 1
                re[it,] = c(idx[i], idx[j])
            }
          }
       }
       colnames(re) = c("id1", "id2")
       return(re)
     }

     re = matrix(NA, ncol=2, nrow=reRow)
     if(re.ordered){
        it = 0
        for ( i in 1:ct ){
          j = i + 1 
          while ( j <= ct){
              ## print(i)
              ## print(j)
              it = it + 1
              re[it,] = range(c(idx[i], idx[j]))
              j = j + 1
          }
        }
     }else{
        it = 0
        for ( i in 1:ct ){
          j = i + 1
          while( j <= ct){
              it = it + 1
              re[it,] = c(idx[i], idx[j])
              j = j + 1
          }
        }
     }
     colnames(re) = c("id1", "id2")
     return(re)

}


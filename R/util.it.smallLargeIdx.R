util.it.smallLargeIdx <-
function(len, keep.same=F){

      if(keep.same){
               init = 1:len
                      idx.ma = matrix(NA, ncol=2 , nrow=(len^2-len)/2+len )
                      idx.init = unlist(lapply(init, FUN=function(i, ct){
                                 rep(i, times=ct-i+1)}, ct=len))
                      idx.next = unlist(lapply(init, FUN=function(i, ct){
                                 i:ct
                               }, ct=len))
                      idx.ma[,1]=idx.init
                      idx.ma[,2]=idx.next

             }else{
                      init=1:(len-1)
                      #print("util.it.smallLargeIdx2")
                      #print(init)
                      #print(len)
                             idx.ma = matrix(NA, ncol=2 , nrow=(len^2-len)/2 )
                             idx.init = unlist(lapply(init, FUN=function(i, ct){
                                        rep(i, times=ct-i)}, ct=len))
                             idx.next = unlist(lapply(init, FUN=function(i, ct){
                                        (i+1):ct
                                      }, ct=len))
                             idx.ma[,1]=idx.init
                             idx.ma[,2]=idx.next
                    }

          return(idx.ma)
    }


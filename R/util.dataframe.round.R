util.dataframe.round <-
function(vecData, keepZero=F, na.replace=NULL, digits=getOption("digits"), ...){

  re = sapply(vecData, FUN=function(x, na.rep, digits, ...){
                                num = x
                                if(is.numeric(num)){

                                    if(is.na(num) & (!is.null(na.rep))){
                                      num = na.replace
                                    }else if(is.na(num) & is.null(na.rep)){
                                      num = NA
                                    }else{
                                      if(keepZero){
                                        if(digits==0){
                                          num = format(num, nsmall=0)
                                        }else{
                                          num = format(num, digits=digits, nsmall=digits)
                                        }
                                      }else{
                                        num = round(num, digits=digits, ...)
                                      }
                                    }
                                                                 
                                }else{
                                  if(num=="NA" & is.null(na.rep)){
                                    num = as.character(num)
                                  }else if(num=="NA" & (!is.null(na.rep))){
                                    num = na.rep
                                  }else{
                                    num = num
                                  }
                                }
                                  
                              }, na.rep = na.replace, digits=digits, ...)
  return(re)

}


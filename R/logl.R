logl <-
function (fileName, str){
       cat(str , file = fileName, sep = " ", fill = T, labels = NULL, append = T)
       return(NULL)
}


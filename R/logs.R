logs <-
function (fileName, str){
       cat(str , file = fileName, sep = " ", fill = F, labels = NULL, append = T)
       return(NULL)
}


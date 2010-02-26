logWarn <-
function (fileName, str){
       cat("************ WARNING ***************** WARNING ***********" , file = fileName, sep = " ", fill = T, labels = NULL, append = T)
       cat(str , file = fileName, sep = " ", fill = T, labels = NULL, append = T)
       cat("************ WARNING ----------------- WARNING ***********" , file = fileName, sep = " ", fill = T, labels = NULL, append = T)
       return(NULL)
}


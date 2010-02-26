logErr <-
function (fileName, str){
       cat("!!!!!!!!!!!! ERROR !!!!!!!!!!!!!!!!!!! ERROR !!!!!!!!!!!!" , file = fileName, sep = " ", fill = T, labels = NULL, append = T)
       cat(str , file = fileName, sep = " ", fill = T, labels = NULL, append = T)
       cat("!!!!!!!!!!!! ERROR ------------------- ERROR !!!!!!!!!!!" , file = fileName, sep = " ", fill = T, labels = NULL, append = T)
       return(NULL)
}


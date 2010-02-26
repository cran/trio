loge <-
function (fileName, str=NULL){
       if(!is.null(str)) cat(str , file = fileName, sep = " ", fill = T, labels = NULL, append = T)
       cat(paste("-------------End of the Script::" , Sys.time(), "--------------------") , file = fileName, sep = " ", fill = T, labels = NULL, append = T)
       return(NULL)
}


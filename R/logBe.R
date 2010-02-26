logBe <-
function (fileName, str=NULL){
       cat(paste("\n\n=============Begin of the Script::" , Sys.time(), "==================") , file = fileName, sep = " ", fill = T, labels = NULL, append = T)
       if(!is.null(str)) cat(str , file = fileName, sep = " ", fill = T, labels = NULL, append = T)
       return(NULL)
}


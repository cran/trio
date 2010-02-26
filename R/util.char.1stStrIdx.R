util.char.1stStrIdx <-
function(str, find){

   result = qstrsplit(str, find, re.1st=T)
   if(length(result)==1){
     # NA==not found it
     if(is.na(result)) return(0)
     # NULL==contain find only
     if(is.null(result)) return(1)
   }

   return(nchar(result[1])+1)
}


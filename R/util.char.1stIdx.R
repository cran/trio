util.char.1stIdx <-
function(str, find, match=T, is.array=F){
   
   if(!is.array){
     charArray = util.str.2CharArray(str, len=nchar(str))
   }else{
     charArray=str
   }
   if(match){
     pos = which(charArray==find)
   }else{
     pos = which(charArray!=find)

     #print(pos)
   }
   if( length(pos)>=1){
     if(match){
       return(pos[1])
     }else{
       return(pos[1]-1)
     }
   }else{
     return (0)
   }
}


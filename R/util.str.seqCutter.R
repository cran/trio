util.str.seqCutter <-
function(str, delims){
   i = 1
   j = 1
   stopS = F
   re = NULL
   len = length(delims)
   curStr = str
   lastIdx = 0
   while (i<=len & !stopS){
     idx.f = util.char.1stStrIdx(curStr, delims[i])
     if (idx.f==0){
       stopS = T
       return(NA)
     }else{
       # find the current one
       #print(curStr)
       #print(idx.f)
       seg = substr(curStr, 1, idx.f-1)
       if(i==len) lastIdx = nchar(curStr)
       curStr = substr(curStr, idx.f+nchar(delims[i]), nchar(curStr))
       #print(curStr)
       re = c(re, seg)
     }
     i = i +1 
   }
   if(idx.f+nchar(delims[i])<lastIdx) re = c(re, curStr) 
   return(re)
 }


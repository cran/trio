util.str.rmSpacePadder <-
function(str, rm=c("b", "e")){

    charArray = util.str.2CharArray(str, len=nchar(str))
    blankPos = which(charArray!=" ")
    if(length(blankPos)==0) return("")
    if(is.element("b", rm)){
       if(is.element("e", rm)){
         re = substr(str, blankPos[1], blankPos[length(blankPos)])
       }else{
         re = substr(str, blankPos[1], nchar(str))
       }
    }else{
       if(is.element("e", rm)){
         re = substr(str, 1, blankPos[length(blankPos)])
       }else{
         re = str
       }
    }
    return(re)
}


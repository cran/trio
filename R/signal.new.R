signal.new <-
function(coef, type, signalStr){

  # var used internally
  dig2Code=c("00", "11", "12", "22")
  
  ifD = F
  vars = NULL
  segReplace = NULL

  
  signal=binaTree.parser(str=signalStr)
  #signal =binaTree.parser.proc(str=signalStr, binaTree=binaTree, curLevel=0, first.seg=T)
  vars = unique(unlist(signal$elm.vlist))
  #boolOp = signal$elm.vMa[,3]

  ori.vars = unlist(signal$elm.vlist)
  ori.match = match(ori.vars, vars)
  
  
  if(type=="geno2d"){
    ## translate into snp 2-digit coding
    pos.equalSign = sapply(vars, FUN=util.char.1stIdx, find="=")

    snpIdx = sapply(1:length(vars), FUN=function(i, name, pos) {
      ## variable name starts with g, followed by snp idx and "="
      a = substr(name[i], 2, pos[i]-1)
      a
    }, name=vars, pos=pos.equalSign)
    
    segReplace = sapply(1:length(vars), FUN=function(i, name, pos, dig2Code, snpId){
      ##  Dec08Change!!! Change the sigStr to two-digit string representation
      a = substr(name[i], pos[i]+1, pos[i]+2)
      ## change!!!change
      if(a== dig2Code[2]){
        ## less common homo is "11", (a) is dominant for the less common, (b) is recessive for the less common
        ## 11 is equivalent to recessive 
        re = paste("v", snpId[i], "b", sep="")  
      }
      if(a== dig2Code[4]){
        ## "22" is equivalent to not dominant
        re =  paste("(not v", snpId[i] ,"a)", sep="")
      }      
      if(a== dig2Code[3]){
        ## heto: "12" 
        re = paste("( (not v", snpId[i],"b) and v", snpId[i], "a )", sep="")    
      }
      re
    }, name=vars, pos=pos.equalSign, dig2Code=dig2Code, snpId = snpIdx)
  } ## if(type=="geno2d"){

  if(type=="D/R"){
    ## translate into snp 2-digit coding
    snpIdx = sapply(1:length(vars), FUN=function(i, name) {
      ## variable name starts with snp idx, followed by one character "D"/"R"
      a = substr(name[i], 1, nchar(name[i])-1)
      a
    }, name=vars)

    snpDR = sapply(1:length(vars), FUN=function(i, name) {
      ## variable name starts with snp idx, followed by one character "D"/"R"
      a = substr(name[i], nchar(name[i]), nchar(name[i]))
      a
    }, name=vars)


    ## need to take care of not
    #if (length(vars)!=length( boolOp )) stop("Error in signal.new: not unique markers in string.")
    
    segReplace = sapply(1:length(vars), FUN=function(i,  snpDR, snpId) {
      
      if(snpDR[i]== "D"){
        ## dominant:
        re = paste( "v", snpId[i] ,"a", sep="")
      }
      if(snpDR[i]== "R"){
        ## recessive:
        re = paste( "v", snpId[i] ,"b", sep="")
      }
      re
    },  snpDR=snpDR,  snpId = snpIdx)
  } ## if(type=="d/r"){
  
    ## replace the original str with new element within str
    changedStr = signalStr
    if(length(vars)>=1){
      for (i in 1:length(vars)){
        tmp.str = vars[i]
        ##print(tmp.str)
        changedStr = util.str.replace(str=changedStr,
                           replaced=vars[i],
                           new=segReplace[i],
                           replace.all=T)
        ##print(model.signal.cur)
      }
    }

  if(ifD) print(changedStr)
  if( length(unlist(signal$elm.vlist))==1 & substr(changedStr,1,1)=="(")
    changedStr = substr(changedStr, 2, nchar(changedStr)-1)
  
  ## reconstruct the single
  #signal =binaTree.parser.proc(str=changedStr, binaTree=binaTree, curLevel=0, first.seg=T)
  signal = binaTree.parser(str=changedStr)
  
  ## also return the snpIdx so reshuffle the column can be done later
  all.varSnp = snpIdx[ ori.match ]
  
  return(list(coef=coef, signal=signal, var.idx = all.varSnp, str=changedStr)) 

}


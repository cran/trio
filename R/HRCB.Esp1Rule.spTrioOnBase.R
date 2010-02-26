HRCB.Esp1Rule.spTrioOnBase <-
function(bkMap=NULL, preObj=NULL, spStrata, rule, caseNo, ifS="simuInfo",  reControl=F){
  FN = "HRCB.Esp1Rule.spTrioOnBase"
  if(is.null(preObj)){
    if(is.null(bkMap)) stop(" No object is passed as bkMap nor as preObj. The function need at least one to work.")
    warning( "No object is passed as preObj. The function will need to generate the preObj.")
    
    preObj = bkMap.HRCB.Esp1Rule.genoSeq(bkMap, rule, re.probOnly = F)
  }

  #print(FN)
  spStrataObj  = get(spStrata)
  #print(str(spStrataObj))

  
  bkMapS=preObj$bkMapS
  supHapExp=preObj$supHapExp
  supHapProb=preObj$supHapProb

  if(is.null(preObj$supHapExp)) stop("Object, preObj, has the variable $supHapExp as NULL!")
  
  kids = HRCB.Esp1Rule.sampleKid(rule=rule, caseNo=caseNo, spStrata=spStrataObj, supHapProb=supHapProb)
    
#  print("return kids")
#  return(kids)
  parIdx = matrix(NA, nrow = caseNo, ncol=4)
  for( eachChild in 1:caseNo){
     hapPair = kids[eachChild,]
     if (hapPair[1]==hapPair[2]){
       parIdx[eachChild,] = ESp.impuParent.homoKid(hapPair[1], supHapProb)
     }else{
       ## testing  hapPair=c(1,2)
       parIdx[eachChild,] = ESp.impuParent.heterKid(hapPair, supHapProb)
     }
   }

    ## shuffle the trio, so risk groups are mixed
    tmpRandomShu = sample(1:caseNo, size=caseNo, replace=F)

   if(!reControl){
     ## need to combine family and run.
     fam.hapIdx.unsort = rbind(parIdx[,1:2], parIdx[,3:4], kids)
     fam.hapIdx = fam.hapIdx.unsort[t( matrix(1:(caseNo*3), ncol=3, nrow=caseNo, byrow=F)), ]
  
     ## find out the exact string expression
     hap.str1 = supHapExp[fam.hapIdx[,1]]
     hap.str2 = supHapExp[fam.hapIdx[,2]]
  
     ## convert string into digits
     geno.FMCMa = covDipStr2CodedGeno(hap.str1, hap.str2, subjectCt=caseNo*3, snpLen=bkMapS$snpLen, snpCoding=0:3, snpBase=c(0, bkMapS$alleleCode))

     ## shuffle the trio, so risk groups are mixed
     allRowIdx = matrix(1:(3*caseNo), nrow=3, byrow=F) 
     trioShuffledIdx = allRowIdx[,tmpRandomShu]
     geno.FMCMa=geno.FMCMa[trioShuffledIdx,]

     if(!is.null(ifS)) {
       ## check!!!
       supDipIdx = matrix(unlist(apply(fam.hapIdx.unsort, 1, FUN=util.it.triMatch, len=length(supHapProb))),
                    ncol=3, nrow=caseNo, byrow=F)
  
       sample.idx = data.frame(parIdx, kids, supDipIdx)[tmpRandomShu,]
       colnames(sample.idx) = c("hapIdx_f1", "hapIdx_f2", "hapIdx_m1", "hapIdx_m2",
                                "hapIdx_c1", "hapIdx_c2", "dipIdx_f", "dipIdx_m", "dipIdx_c")
       write.csv(sample.idx,  file=paste(ifS, "supHap.csv", sep=""))
     }
   }else{ # if(!reControl){
     ## need to find out which one is for the affected child
     othKids =  matrix(NA, nrow = caseNo, ncol=6)
     for( eachChild in 1:caseNo){
       find.idx1 = which( kids[eachChild, 1] == parIdx[eachChild, 1:4])
       if(length(find.idx1)==0) stop("Affected child hap indexes do not match with parents'.")
       find.idx2 = which( kids[eachChild, 2] == parIdx[eachChild, 1:4])
       if(length(find.idx2)==0) stop("Affected child hap indexes do not match with parents'.")
   
       if(max(find.idx1)<=2){
         # idx1 can only come from F, then idx2 must come from M
         fix.idx1 = find.idx1[1]
         fix.idx2 = find.idx2[ find.idx2>=3 ][1]
         comp.idx1 = (1:2)[-fix.idx1]
         comp.idx2 = (3:4)[-(fix.idx2-2)]
       }else{ ## if(max(find.idx1)<=2){
         if( min(find.idx1)>=3){
           # idx1 can only come from M, then idx2 must come from F
           fix.idx1 = find.idx1[1]
           fix.idx2 = find.idx2[ find.idx2<=2 ][1]
           comp.idx1 = (3:4)[-(fix.idx1-2)]
           comp.idx2 = (1:2)[-fix.idx2]
         }else{ ### if( min(find.idx1)>=3){
           # idx1 can come from either F or M, then need to see where idx2 come from
           if(max(find.idx2)<=2){
             # idx2 must come from F, then idx1 come from M
             fix.idx1 = find.idx1[ find.idx1>=3 ][1]
             fix.idx2 = find.idx2[1]
             comp.idx1 = (3:4)[-(fix.idx1-2)]
             comp.idx2 = (1:2)[-fix.idx2] 
           }else{ #### if(max(find.idx2)<=2){
#             if(min(find.idx2)>=3){
#               # idx2 must come from M, then idx1 come from F
#               fix.idx1 = find.idx1[ find.idx1<=2 ][1]
#               fix.idx2 = find.idx2[1]
#               comp.idx1 = (1:2)[-fix.idx1]
#               comp.idx2 = (3:4)[-fix.idx2]
#             }else{
               # set idx1 come from F, then idx2 come from M
               fix.idx1 = find.idx1[ find.idx1<=2 ][1]
               fix.idx2 = find.idx2[ find.idx2>=3 ][1]
               comp.idx1 = (1:2)[-fix.idx1]
               comp.idx2 = (3:4)[-(fix.idx2-2)]               
#             }
           } #### if(max(find.idx2)<=2){
         } ### if( min(find.idx1)>=3){
       } ## if(max(find.idx1)<=2){
       # fill in the other three
       
       if(  max( is.na(c(fix.idx1, comp.idx2, comp.idx1, comp.idx2, comp.idx1, fix.idx2)))==1 ) print("##############")
       othKids[eachChild, ] =  parIdx[eachChild,1:4] [c(fix.idx1, comp.idx2, comp.idx1, comp.idx2, comp.idx1, fix.idx2)]
       #print( c(parIdx[eachChild,1:4], kids[eachChild,], othKids[eachChild,]) )
       
     } # for( eachChild in 1:caseNo){

     ## need to combine family and run.
     fam.hapIdx.unsort = rbind(parIdx[,1:2], parIdx[,3:4], kids, othKids[,1:2], othKids[,3:4], othKids[,5:6])
     fam.hapIdx = fam.hapIdx.unsort[t( matrix(1:(caseNo*6), ncol=6, nrow=caseNo, byrow=F)), ]
  
     ## find out the exact string expression
     hap.str1 = supHapExp[fam.hapIdx[,1]]
     hap.str2 = supHapExp[fam.hapIdx[,2]]
  
     ## convert string into digits
     geno.FMCMa = covDipStr2CodedGeno(hap.str1, hap.str2, subjectCt=caseNo*6, snpLen=bkMapS$snpLen, snpCoding=0:3, snpBase=c(0, bkMapS$alleleCode))

     ## shuffle the trio, so risk groups are mixed
     allRowIdx = matrix(1:(6*caseNo), nrow=6, byrow=F) 
     trioShuffledIdx = allRowIdx[,tmpRandomShu]
     geno.FMCMa=geno.FMCMa[trioShuffledIdx,]

     if(!is.null(ifS)) {
       ## check!!!
       supDipIdx = matrix(unlist(apply(fam.hapIdx.unsort, 1, FUN=util.it.triMatch, len=length(supHapProb))),
                    ncol=6, nrow=caseNo, byrow=F)
  
       sample.idx = data.frame(parIdx, kids, supDipIdx)[tmpRandomShu,]
       colnames(sample.idx) = c("hapIdx_f1", "hapIdx_f2", "hapIdx_m1", "hapIdx_m2",
                                "hapIdx_c1", "hapIdx_c2", "dipIdx_f", "dipIdx_m", "dipIdx_c", "dipIdx_cA", "dipIdx_cB" ,"dipIdx_cC" )
       write.csv(sample.idx,  file=paste(ifS, "supHap.csv", sep=""))
     }
   } # if(!reControl){
  
     
   return(geno.FMCMa)

}


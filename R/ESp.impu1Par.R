ESp.impu1Par <-
function(othParPairs, childPairs,  semiMapFrame, resiProbCol, augIdxCol, probCol, snpLen, reType=F, job=1 ){
  ifD = F
  fStr = "[ESp.impu1Par]"
  if(ifD) {
    print(fStr)
    print(semiMapFrame)
    print(paste("snpLen=", snpLen))
  }
  
  if(is.null(dim(othParPairs))) othParPairs = matrix(othParPairs, ncol=2)
  if(is.null(dim(childPairs))) childPairs = matrix(childPairs, ncol=2)
  
  allParIdx= unique(as.vector(othParPairs))

  if(ifD){
    print(othParPairs)
    print(childPairs)

  }
  
  
  ## generate all combination of par and child pairs.
  prob.ps = getHapProb2(selIdx=allParIdx,
    semiMapFrame,
    resiProbCol=resiProbCol,
    augIdxCol=augIdxCol,
    probCol= probCol,
    snpLen, restandard=T)
  if(ifD) {
    print(prob.ps)
  }
  
  n.par = length(othParPairs)/2
  n.ch = length(childPairs)/2

  sp.ma = NULL
  for( i in 1:n.par){
    par = othParPairs[i,]
    
    for( j in 1:n.ch){
      ch = childPairs[j,]
      
      ## find whether the pair match on at lease one hap
      if(sum(is.na(match(par, ch)))==2) {
        if(ifD) print( paste("Ignor par=", paste(par, collapse=";"), " and ch=", paste(ch, collapse=";") ))

      }else{
          sp.prob=NULL
          type = 0
          
          ##1) B/E/C vs. A/D
          if (par[1]!=par[2]){
            ##2) B/E/C, B/E vs. C
            if(ch[1]!=ch[2]){
              
              ch = sort(ch)
              par= sort(par)
              ##3) B/E, B vs. E
              if (sum( ch==par )==2){
                if(ifD) print("B")
                # B
                sp.prob = ESp.impu1Par.B(hap=ch, prob.p=prob.ps[match(par, allParIdx)], semiMapFrame, resiProbCol, augIdxCol, probCol, snpLen )
                type=2

                ## additional col, no use
                add.col = rep(NA, 5)
                
              }else{
                if(ifD) print("E")
                # E: need to keep the common one at the front
                match.t = is.na(match(ch, par))
                if(match.t[1]){ ch = c(ch[2], ch[1])}
                match.t = is.na(match(par, ch))
                if(match.t[1]){ par = c(par[2], par[1]) }             
                
                sp.prob = ESp.impu1Par.E(hap=ch[1], prob.p=prob.ps[match(par, allParIdx)], semiMapFrame, resiProbCol, augIdxCol, probCol, snpLen )
                type=5
                add.col = rep(par[2], 2)
                if(ifD) print(sp.prob)
                
              } ##B vs E::if (sum( ch==par )==2){
              
            }else{ ##2) B/E/C, B/E vs. C
              # C: need to keep the common one at front
              if(ifD) print("C")
              match.t = is.na(match(par, ch))
              if(match.t[1]){par = c(par[2], par[1]) }

              #print(par)
              #print( prob.ps[match(par, allParIdx)] )

              sp.prob = ESp.impu1Par.C(hap=par, prob.p=prob.ps[match(par, allParIdx)], semiMapFrame, resiProbCol, augIdxCol, probCol, snpLen )
              # print(sp.prob)
              type=3
              add.col = rep(par[2], 3)

              
            } ## B/E vs C::if(ch[1]!=ch[2]){
            
          }else{ ## B/E/C vs D/A if (par[1]!=par[2]){
            ##2) D/A, D vs. A
            if(ch[1]!=ch[2]){
              # D need to keep the common one at front
              if(ifD) print("D")
              match.t = is.na(match(ch, par))
              if(match.t[1]){
                  ch = c(ch[2], ch[1])              
              }
              sp.prob = ESp.impu1Par.D(hap=ch[2], prob.p=prob.ps[match(par[1], allParIdx)], semiMapFrame, resiProbCol, augIdxCol, probCol, snpLen )
              type=4
              
              add.col = rep(NA, 2)
            }else{
              # A
              if(ifD) print("A")
              sp.prob = ESp.impu1Par.A(hap=ch[1], prob.p=prob.ps[match(par[1], allParIdx)], semiMapFrame, resiProbCol, augIdxCol, probCol, snpLen )
              type=1
              add.col = rep(NA, 2)
              
            } ## D vs. A::if(ch[1]!=ch[2]){
    
          } ## B/E/C vs D/A if (par[1]!=par[2]){
          ch.colA = matrix(rep(ch, length(sp.prob)), ncol=2, byrow=T)
          
          rows = cbind(rep(type, length(sp.prob)), 1:length(sp.prob), sp.prob, add.col, ch.colA)
          if(ifD) print(rows)
          
          sp.ma = rbind(sp.ma, rows)
          # print(sp.ma)

        } # if(sum(is.na(match(par, ch)))<2) {
    } # for( j in 1:n.ch){
  } # for( i in 1:n.par){

  colnames(sp.ma)=c("type", "straSeq", "prob", "add", "ch1", "ch2")


  hap6idx = matrix(NA, nrow=job, ncol=7)
  
  for(ss in 1:job){
    # now sample it
    row.sp = sample(1:nrow(sp.ma), size=1, prob=sp.ma[,3])
  
    type.sp = sp.ma[row.sp,1]
    straSeq = sp.ma[row.sp,2]
    
    if(ifD) {
      print("Sampling")
      print(sp.ma)
    }
    if(type.sp==1){
      sp6 = ESp.impu1Par.A.sp(hap=sp.ma[row.sp, 5],  straSeq, semiMapFrame, resiProbCol, augIdxCol, probCol, snpLen )
    }
    if(type.sp==2){
      sp6 = ESp.impu1Par.B.sp(hap=sp.ma[row.sp, 5:6],  straSeq, semiMapFrame, resiProbCol, augIdxCol, probCol, snpLen )
    }
    if(type.sp==3){
      #print("type.sp==3")
      sp6 = ESp.impu1Par.C.sp(hap=sp.ma[row.sp, 5], hap.p=sp.ma[row.sp, 4], straSeq, semiMapFrame, resiProbCol, augIdxCol, probCol, snpLen )
      #print(sp6)
    }
    if(type.sp==4){
      sp6 = ESp.impu1Par.D.sp(hap=sp.ma[row.sp, 5:6],  straSeq, semiMapFrame, resiProbCol, augIdxCol, probCol, snpLen )
    }
    if(type.sp==5){
      sp6 = ESp.impu1Par.E.sp(hap=sp.ma[row.sp, 5:6], hap.p=sp.ma[row.sp, 4] ,  straSeq, semiMapFrame, resiProbCol, augIdxCol, probCol, snpLen )
    }

    hap6idx[ss, ]=c(sp6, type.sp)
  }

  
  if(reType){
    return(hap6idx )
  }else{
    return(hap6idx[,1:6,drop=F])
  }
}


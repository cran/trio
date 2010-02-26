binaTree.toStr <-
function(binaTree){
  
  ifD = F

  ## process the variable only element first, regardless the level
  vCt = nrow(binaTree$elm.vMa)
  data.v = rep("", times=vCt)
  
  for( i in 1:vCt){
    # print(i)
    cur.row = binaTree$elm.vMa[i,]
    ## col seq= id, level, bool
    ## find the var list
    cur.elm = binaTree$elm.vlist[[ i ]]

    cur.bool = cur.row[3]
    if(length(cur.elm)>1){
     
      if(cur.bool==0){ ## and
        data.v[i]= paste(cur.elm, collapse=" and ")
      }
      if(cur.bool==1){ ## or
        data.v[i]= paste(cur.elm, collapse=" or ")
      }
    }else{
      if(cur.bool==-1){ ## not
        data.v[i]= paste("(not ", cur.elm, ")", sep="")
      }else{
        data.v[i]= cur.elm
      }
    }
  } # for( i in 1:vCt){

  ## if have node element, process the node from deeper level to zero
  data.n=NULL
  nCt = nrow(binaTree$elm.nMa)
  if(nCt>=1){
    nMa = binaTree$elm.nMa
    
    data.n=rep("", times=nCt)

    new.order = order(nMa[,2], decreasing = TRUE)
    if(ifD) print(nMa)
    nMa = nMa[new.order, ,drop=F]
    if(ifD) print(nMa)
    #nlist = binaTree$elm.nlist[new.order]
    nlist = binaTree$elm.nlist
    if(ifD) print(nlist)
    for( i in 1:nCt ){
      cur.row = nMa[i, ,drop=F]
      if(ifD) print("cur node row after sorting")
      if(ifD) print(cur.row)
      ## col seq= id, level, bool, open
      
      ## find the var list
      cur.elm = nlist[[ cur.row[1] ]]
      if(ifD) print("cur node element")
      if(ifD) print(cur.elm)
      dd.velm=NULL
      dd.nelm=NULL
      # if variable element, get from data.v
      if(sum(cur.elm[,2]==0)>=1){
        dd.velm = data.v[cur.elm[,1][cur.elm[,2]==0]]
      }
      # if node element, get from data.n
      if(sum(cur.elm[,2]==1)>=1){
        dd.nelm = data.n[cur.elm[,1][cur.elm[,2]==1]]
      }

      
      dd.cur = cbind(dd.velm, dd.nelm)
      if(ifD) print(cur.row)
      cur.bool = cur.row[3]

      if(i!=nCt){
        if(cur.bool==0){ ## and
          data.n[cur.row[1]]= paste("(" , paste(dd.cur, collapse= " and "), ")", sep="")
        }
        if(cur.bool==1){ ## or
          data.n[cur.row[1]]= paste("(" , paste(dd.cur, collapse= " or "), ")", sep="")
        }
      }else{ # for the top node, no need for outside brace
        if(cur.bool==0){ ## and
          data.n[cur.row[1]]= paste(dd.cur, collapse= " and ")
        }
        if(cur.bool==1){ ## or
          data.n[cur.row[1]]= paste(dd.cur, collapse= " or ")
        }  
      }
    } ## for( i in 1:nCt ){
    
  } ##  if(nCt>=1){

    if(nCt==0){
      # only one variable only element
      return(data.v[,1])
    }else{
      # return the last element (must be a node)
      t.id = binaTree$curElmId
      re = data.n[t.id]
      if(ifD) print(data.n)
      return(re)
    }

}


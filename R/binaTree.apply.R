binaTree.apply <-
function(binaTree, data, colnames, re.final=T){
  ifD = F
  ori.rowCt = nrow(data)
  col.idx = 1:ncol(data)
  
  ## data.v: a matrix of boolean (0/1) for each row in binaTree$vMa, ordered by id
  ## data.n: a matrix of boolean (0/1) for each row in binaTree$nMa, ordered by id
  data.v=NULL
  data.n=NULL

  ## data.f: a vector of boolean (0/1) for the node at level 0
  ## or for the variable element if no node is defined.
  data.f=NULL

  data.f = rep(NA, ori.rowCt)
  
  ## process the variable only element first, regardless the level
  vCt = nrow(binaTree$elm.vMa)
  data.v = matrix(NA, ncol=vCt, nrow=ori.rowCt) 
  for( i in 1:vCt){
    cur.row = binaTree$elm.vMa[i,]
    ## col seq= id, level, bool
    ## find the var list
    cur.elm = binaTree$elm.vlist[[ i ]]
    cur.col = col.idx[is.element(colnames, cur.elm)]
    if (length(cur.col)!=length(cur.elm)) stop()

    cur.bool = cur.row[3]
    if(length(cur.col)>1){
     
      if(cur.bool==0){ ## and
        data.v[,i]=as.numeric(apply(data[,cur.col], 1, all))
      }
      if(cur.bool==1){ ## or
        data.v[,i]=as.numeric(apply(data[,cur.col], 1, any))
      }
    }else{
      if(cur.bool==-1){ ## not
        data.v[,i]=as.numeric(1-data[,cur.col])
      }else{
        data.v[,i]=data[,cur.col]
      }
    }
  } # for( i in 1:vCt){

  ## if have node element, process the node from deeper level to zero
  data.n=NULL
  nCt = length(binaTree$elm.nlist)
  if(nCt>=1){
    nMa = binaTree$elm.nMa
    
    data.n=matrix(NA, ncol=nCt, nrow=ori.rowCt)

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
        dd.velm = data.v[, cur.elm[,1][cur.elm[,2]==0]]
        
        if (ifD) {
          if(is.null(dim(dd.velm))){
            print(dd.velm[1:5])
          }else{
            print(dd.velm[1:5,])
          }
        }
      }
      # if node element, get from data.n
      if(sum(cur.elm[,2]==1)>=1){
        dd.nelm = data.n[, cur.elm[,1][cur.elm[,2]==1]]
        
        if (ifD) {
          if(is.null(dim(dd.nelm))){
            print(dd.nelm[1:5])
          }else{
            print(dd.nelm[1:5,])
          }
        }
      }

      dd.cur = cbind(dd.velm, dd.nelm)
      if(ifD) print(cur.row)
      cur.bool = cur.row[3]
     
      if(cur.bool==0){ ## and
        data.n[,cur.row[1]]=as.numeric(apply(dd.cur, 1, all))
      }
      if(cur.bool==1){ ## or
        data.n[,cur.row[1]]=as.numeric(apply(dd.cur, 1, any))
      }

      if( (cur.bool!=1) & (cur.bool!=0)) stop("Wrong boolean value.")
    } ## for( i in 1:nCt ){
    
  } ##  if(nCt>=1){

  if(re.final){
    if(nCt==0){
      # only one variable only element
      return(data.v[,1])
    }else{
      # return the last element (must be a node)
      t.id = binaTree$curElmId
      re = data.n[,t.id]
      return(re)
    }

  }else{
    return(list(data.v=data.v, data.n=data.n))
  }
}


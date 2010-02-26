binaTree.setApply <-
function(binaTree, setList){
  ifD = F
  # the setList already contain the set for each variable

  data.n=NULL
  
  nCt = length(binaTree$elm.nlist)
  if(nCt>=1){
    nMa = binaTree$elm.nMa
    
    data.n = rep(list(NULL), nrow(nMa))

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
      ## only three possibilties: both are variable, both are node, one variable and one node
      if( sum(cur.elm[,2]==0)  >=1){
        # dd.velm = data.v[, cur.elm[,1][cur.elm[,2]==0]]
        dd.velm = setList[ c(cur.elm[,1][cur.elm[,2]==0]) ]
        
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
        dd.nelm = data.n[ c(cur.elm[,1][cur.elm[,2]==1]) ]
        
        if (ifD) {
          if(is.null(dim(dd.nelm))){
            print(dd.nelm[1:5])
          }else{
            print(dd.nelm[1:5,])
          }
        }
      }

      dd.cur = c(dd.velm, dd.nelm)
      if(ifD) print(cur.row)
      cur.bool = cur.row[3]
     
      if(cur.bool==0){ ## and -> intersection

        if(length( dd.cur[[1]] )==0) {
          warning("length( dd.cur[[1]] )==0")
          data.n[[ cur.row[1]   ]] = NULL
        }
        if(length( dd.cur[[2]] )==0) {
          warning("length( dd.cur[[2]] )==0")
          data.n[[ cur.row[1]   ]] = NULL
        }
        if(length(dd.cur[[1]])!=0 & length(dd.cur[[2]])!=0){
          filter = match(dd.cur[[1]], dd.cur[[2]], nomatch=0)
          filter = filter[filter>0]
          data.n[[  cur.row[1]  ]] = dd.cur[[2]][ filter ]
        }
      }
      if(cur.bool==1){ ## or -> union
        
        data.n [[  cur.row[1]  ]] = unique(c(dd.cur[[1]], dd.cur[[2]]))
        if( length(dd.cur[[1]])==0 & length(dd.cur[[2]])==0 ){
          data.n[[  cur.row[1]  ]] = NULL
        }
      }
    } ## for( i in 1:nCt ){
    
  } ##  if(nCt>=1){

    if(nCt==0){
      # only one variable only element
      return(setList[[1]])
    }else{
      # return the last element (must be a node)
      t.id = binaTree$curElmId
      re = data.n[[t.id]]
      return(re)
    }

}


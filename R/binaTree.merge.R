binaTree.merge <-
function(treeForSearch){
 fN = "binaTree.merge:"
 ifD = F
 
 searchWholeTree = TRUE
 tt = 0
 search = FALSE
 # do an iterative search for ( or ) and , keep search until all the "or" have no parent with "and" operator
 while(searchWholeTree){
     tt = tt+1
     if(tt>100) stop(paste(fN, "iterative procdure exceed 100 times. Sometime is wrong!"))
     
     # if the tree is updated, use the updated one. Otherwise, search the next possible choice
     if(!search){
       if(ifD) {
         print("Updated the tree:")
         #print(binaTree.toStr(treeForSearch))
         print(treeForSearch)
       }
       tree = treeForSearch
       elm.nMa = tree$elm.nMa
  
       ## search for "and"
       filter = elm.nMa[,3]==0
       if(sum(filter)==0) {
  #          search = F
  #          searchWholeTree = F
           return(treeForSearch)
       }
       search = TRUE
     }

     if( sum(filter)>=1){
       
       newMa = elm.nMa[filter, , drop=FALSE]

       if(ifD) print(newMa)
       
       j = 1 # loop though node with "and" operator
       remove.list = NULL
       while ( j <= sum(filter) & search){
         if(ifD) print(paste("j=", j))
         curCh.idx = newMa[j, 1]
         curChRow = match(curCh.idx, elm.nMa[,1])
         curNode.level = newMa[j, 2]
         par.idx = binaTree.shuffle.findNodePar(tree, curCh.idx)

         if(is.null(par.idx)) stop (paste(fN, " No parent for a child node!"))

         if(par.idx != 0) {

            parRow = match(par.idx, elm.nMa[,1])
            if(ifD){
              
              print(paste("Matched parent row =", parRow, " parent id =", par.idx))
            }
            # find parent node is using "and"
            if( elm.nMa[parRow, 3] == 0){
              if(ifD) print("find a need for merge")
              search = FALSE
       
              ## find out the element of current node and the other node id
              curNodeMa = tree$elm.nlist[[curCh.idx]]
       
              parNodeMa = tree$elm.nlist[[par.idx]]

              ## need to be the same id as well a node
              filter.curkidRow = parNodeMa[,1]==curCh.idx & parNodeMa[,2]==1

              othNode.row = parNodeMa[!filter.curkidRow,]
              curNode.rows = curNodeMa
              tree$elm.nlist[[par.idx]] = rbind(othNode.row, curNode.rows)

              ## remove from the nMa the current node, will NOT remove the node from nlist for simplicity
              remove.list = c(remove.list, curChRow)

              if(ifD) {
                print("elm.nlist and nMa after updating the partent of old node")
                #print(tree$elm.nlist)
                #print(tree$elm.nMa)
              }
   
              treeForSearch = tree
              tree = NULL
            
            } ## if( elm.nMa[match.id, 3] == 0){
          } # if(par != 0) {
         j = j +1
       } ##  while ( j <= sum(filter) & search){

       if(length(remove.list)>0){
         elm.nMa = elm.nMa[-c(remove.list),,drop=F]
         print(paste("removed:", paste(c(remove.list), collapse=";")))
         if(ifD) print(elm.nMa)
         treeForSearch$elm.nMa = elm.nMa
# 
#          remove.nlist  = elm.nMa[remove.list,1]
#          treeForSearch$elm.nlist = treeForSearch$elm.nlist [-c(remove.nlist)]
       }

       if(j > sum(filter) & search) searchWholeTree = FALSE
       if(ifD & j > sum(filter) & search ) print("stopped")
     } ## if( sum(filter)>=1){
     
   } ## while(searchWholeTree){

 ##return(NULL)
 return(treeForSearch)
}


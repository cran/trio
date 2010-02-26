binaTree.shuffleNode <-
function(treeForSearch){
 fN = "binaTree.shuffleNode:"
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
         print(treeForSearch)
       }
       tree = treeForSearch
       elm.nMa = tree$elm.nMa
  
       ## search for "or"
       filter = elm.nMa[,3]==1
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
       
       j = 1 # loop though node with "or" operator
       while ( j <= sum(filter) & search){
         if(ifD) print(paste("j=", j))
         curCh.idx = newMa[j, 1]
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
              if(ifD) print("find a case")
              search = FALSE
       
              ## find out the element of current node and the other node id
              curNodeMa = tree$elm.nlist[[curCh.idx]]
       
              parNodeMa = tree$elm.nlist[[par.idx]]

              filter.curKid = parNodeMa[,1]==curCh.idx & parNodeMa[,2]==1
              
              twoKid = parNodeMa[,1]
              othNode.idx = twoKid[!filter.curKid]
              othNode.type = parNodeMa[!filter.curKid, 2]

              ## if the other node is a nodel
              ## change the other kid level one down
              if(othNode.type==1){
                othNodeRow = match(othNode.idx, elm.nMa[,1])
                elm.nMa[othNodeRow, 2] =  curNode.level+1

                ## create a duplicate of the othNode
                othNodeMaDup  = tree$elm.nlist[[othNode.idx]]
                othNodeDup.New = length(tree$elm.nlist)+1
                tree$elm.nlist = c(tree$elm.nlist, list(othNodeMaDup))
                elm.nMa = rbind(elm.nMa, c(othNodeDup.New, curNode.level+1, elm.nMa[othNode.idx, 3], 0))
                ma =  matrix(c(curNodeMa[2,1], othNodeDup.New, curNodeMa[2,2], othNode.type), ncol=2, nrow=2, byrow=F)
              }else{
                ma =  matrix(c(curNodeMa[2,1], othNode.idx, curNodeMa[2,2], othNode.type), ncol=2, nrow=2, byrow=F)
              }
              
              ## create the new node
              othNode.New = length(tree$elm.nlist)+1
              tree$elm.nlist = c(tree$elm.nlist, list(ma))
   
              if(ifD) {
                print("elm.nlist after adding the new node")
                #print(tree$elm.nlist)
   
              }
   
              ## added the new node to the nMa
              elm.nMa = rbind(elm.nMa, c(othNode.New, curNode.level, 0, 0))
              
       
              ## update the old node with "and"
              curNodeMa[2,] = c(othNode.idx, othNode.type)
              elm.nMa[curCh.idx, 3]=0
              tree$elm.nlist[[curCh.idx]]=curNodeMa
   
              if(ifD){
                print("elm.nlist after updating the current node")
                #print(tree$elm.nlist)
              }
              
       
              ## update the parent of old node
              parNodeMa =  matrix(c(curCh.idx, othNode.New, 1, 1), ncol=2, nrow=2, byrow=F)
              elm.nMa[parRow, 3]=1
              tree$elm.nlist[[par.idx]]=parNodeMa
              tree$elm.nMa = elm.nMa
              
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
       if(j > sum(filter) & search) searchWholeTree = FALSE
       if(ifD & j > sum(filter) & search ) print("stopped")
     } ## if( sum(filter)>=1){
     
   } ## while(searchWholeTree){

 ##return(NULL)
 return(treeForSearch)

}


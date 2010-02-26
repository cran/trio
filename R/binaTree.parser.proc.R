binaTree.parser.proc <-
function(str, binaTree, curLevel, first.seg=F){
  ## need to grow the list
  ifD = F
  
  str = util.str.rmSpacePadder(str, rm=c("b", "e"))
  
  if(ifD) print(qp("binaTree.parser.proc::str=", str, "|end"))
  if(ifD) print(qp("binaTree.parser.proc::curLevel=", curLevel, "|end"))
  #if(ifD) print(binaTree)
  
   ## patterns: parse out the seg separated by (/)
   ## "" hides the part has already been processed, ')' hides the part not necessary there
   ## 7 "("A and B and C')'
   ## 5  and A and B')'
   ## 8 ")") and A 
   ## 1 "((A and "(B and C)) and
   ## 2   A and (
   ## 4  and (
   
   b.level = util.findBrace(str)
   if(ifD) print(b.level)
   ## if no more brace, it is ending!!!
   t.elm=NULL
   if(is.null(b.level)){

     if(first.seg){
       ## #7:"("A and B and C')' only one elm
       ## can only be and/or/not
       ## add elm, elment only
       if(ifD) print("first.seg")
       t.elm = binaTree.parser.pVarElm(str, type="middle")
       if(is.null(t.elm)) stop()

       ## changed, so it would not allow two elments as one element variable
       if(t.elm$mixed==2){
         if(ifD) print("t.elm$mixed==2")
         # add the node
         
         ## add the element and add the node as well. Note that at the maximum, only two elements is allowed
         binaTree = binaTree.addElm(binaTree, curLevel, t.elm$idx[1], 0)
         binaTree = binaTree.addNode(binaTree, curLevel, t.elm$bool)
         binaTree = binaTree.addElm(binaTree, curLevel, t.elm$idx[2], 0)
         ## close the node
         binaTree = binaTree.closeNode(binaTree, curLevel)
       }else{       
         binaTree = binaTree.addElm(binaTree, curLevel, t.elm$idx, t.elm$bool)
       }
       
       return(binaTree)
     }else{
       ## #5: and A and C')', ending
       ## can be and/or,  must be mix node and
       ## ENDING
       ## add elm, add node, close node
       t.elm = binaTree.parser.pVarElm(str, type="start")
       if(is.null(t.elm)) stop("binaTree.parser.pVarElm(str, type='start')")
       if(t.elm$mixed==1){
         ## add the node
         binaTree = binaTree.addNode(binaTree, curLevel, t.elm$bool)
         ## add elm first
         binaTree = binaTree.addElm(binaTree, curLevel, t.elm$idx, t.elm$bool)
         ## close the node
         binaTree = binaTree.closeNode(binaTree, curLevel)
       }else{
         stop("binaTree.parser.pVarElm(str, type='start'):: not node" )
       }
       return(binaTree)
     } # if(first.seg){

   } # if(is.null(b.level)){
  
   ## if more brace, keep subsetting
   if(!is.null(b.level)){
     if (b.level$brace==1){
       ## #8 or #5 or #7 or having ")"
       if(first.seg){
         stop(" Can not have ')' before '('.")
       }
       
       if(b.level$idx==1){
           ## #8
           ## process node, then subset
           ## close the node
           
           binaTree = binaTree.closeNode(binaTree, curLevel)
           curLevel = curLevel-1
           
           ## subset
           if(b.level$id < nchar(str)){
             curStr = substr(str, b.level$idx+1, nchar(str))
             binaTree = binaTree.parser.proc(curStr, binaTree, curLevel, first.seg=F)
           }
           return(binaTree)
       }else{ # if(b.level$idx==1){
         ## #7 or #5
         ## #7: "("A and B and C')', #5:  and A and B')'
         ## depends on the first is boolean or factor
         ## process node, then subset
         segStr = substr(str, 1, b.level$idx-1)
         t.elm = binaTree.parser.pVarElm(segStr, type="un")
         if(ifD) print(t.elm) 
         if(is.null(t.elm)) stop( "binaTree.parser.pVarElm(str, type='un')" )
         if(t.elm$mixed==1){
           ## add the node
           binaTree = binaTree.addNode(binaTree, curLevel, t.elm$bool)           
           ## print(t.elm)
           ## add elm first
           binaTree = binaTree.addElm(binaTree, curLevel, t.elm$idx, t.elm$bool)
           
           binaTree = binaTree.closeNode(binaTree, curLevel)
           curLevel = curLevel-1
         }else{
           if(t.elm$mixed==2){
             if(ifD) print("second:: t.elm$mixed==2")
             # add the node
             
             ## add the element and add the node as well. Note that at the maximum, only two elements is allowed
             binaTree = binaTree.addElm(binaTree, curLevel, t.elm$idx[1], 0)
             binaTree = binaTree.addNode(binaTree, curLevel, t.elm$bool)
             binaTree = binaTree.addElm(binaTree, curLevel, t.elm$idx[2], 0)
             ## close the node
             binaTree = binaTree.closeNode(binaTree, curLevel)
             curLevel = curLevel - 1
           }else{
               
             binaTree = binaTree.addElm(binaTree, curLevel, t.elm$idx, t.elm$bool)
             ##print(" binaTree.addElm(binaTree, curLevel, t.elm$idx, t.elm$bool) "  )
             ##print(binaTree)
             curLevel = curLevel-1
             ## if the last is an variable only element, need to close the node
           }
         }

         ## subset
         if(b.level$id < nchar(str)){
             curStr = substr(str, b.level$idx+1, nchar(str))
             binaTree = binaTree.parser.proc(curStr, binaTree, curLevel, first.seg=F)
         }else{
             if( curLevel==0){
               binaTree = binaTree.closeNode(binaTree, curLevel)
             }else{
               stop(qp("curLevel=", curLevel, ", should be 0"))
             }
         }
         return(binaTree)
       } # if(b.level$idx==1){
     } # if (b.level$brace==1){

     if(b.level$brace==0){
       ## having "("
       if(b.level$idx==1){
         ## #1: "((A and "(B and C)) and
         ## see how many levels it opens, then subset
           add.lev = util.char.1stIdx(str, "(", match=F)
           curLevel = curLevel+add.lev
           curStr = substr(str, add.lev+1, nchar(str))
           binaTree = binaTree.parser.proc(curStr, binaTree, curLevel, first.seg=F)
           return(binaTree)
       }else{ # if(b.level$idx==1){
         ## #4: and (, #2: A and (

         ## process the node, then subset
         segStr = substr(str, 1, b.level$idx-1)
         t.elm = binaTree.parser.pVarElm(segStr, type="end")
         ##print(paste("t.elm=", t.elm))
         if(is.null(t.elm)) stop( " binaTree.parser.pVarElm(str, type='end')" )
         if(t.elm$mixed==1){
           ## if the boolean term itself
           if(is.null(t.elm$idx)){
             binaTree = binaTree.addNode(binaTree, curLevel, t.elm$bool)
           }else{
             binaTree = binaTree.addElm(binaTree, curLevel, t.elm$idx, t.elm$bool)
             binaTree = binaTree.addNode(binaTree, curLevel, t.elm$bool)
           }

         }else{
           stop( " binaTree.parser.pVarElm(str, type='end'):: not node")        
         }           
         
         ## subset
         if(b.level$id < nchar(str)){
             curStr = substr(str, b.level$idx, nchar(str))
             binaTree = binaTree.parser.proc(curStr, binaTree, curLevel, first.seg=F)
         }else{
           stop ("find ( without )")
         }
         return(binaTree)
       }  # if(b.level$idx==1){
     }# if(b.level$brace==0){
   } # if(!is.null(b.level)){

}


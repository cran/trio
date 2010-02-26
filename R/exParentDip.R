exParentDip <-
function(childTb, parTb, pHapBkInfoMap, snpLen){
	## if a haplotype doesn't include in the parents or child, exclude
	
	if(is.null(dim(childTb)) | is.null(dim(parTb))  | length(pHapBkInfoMap$missingIn)==snpLen){
		return(list(childTb=matrix(childTb, ncol=2), parTb=matrix(parTb, ncol=2)))
	}
	
	child = unique(childTb)
	parent = unique(parTb)
	
	childLeft = child[is.element(child, parent)]
	parentLeft = parent[is.element(parent, child)]
	
	if(length(childLeft)==0 | length(parentLeft)==0 ) {
		## print(paste("ChildTb:\n", childTb, "parTb:\n", parTb, "\n Cannot find the pair matched"))
		stop(paste("\n ChildTb:\n", childTb, "parTb:\n", parTb, "\n Cannot find the pair matched"))
	}
	
	
	childIdx = apply(childTb, 1, FUN = function(rowItem, bench){
				re = as.logical(max(is.element(rowItem, bench)))
			}, bench = childLeft)
	
	
	parIdx = apply(parTb, 1, FUN = function(rowItem, bench){
				re = as.logical(max(is.element(rowItem, bench)))
			}, bench = parentLeft)
	
	if(sum(childIdx)==0 | sum(parIdx)==0 ) {
		## print(paste("childTb:\n", childTb, "parTb:\n", parTb, "\n Cannot find the pair matched"))
		stop(paste("\n childTb:\n", childTb, "parTb:\n", parTb, "\n Cannot find the pair matched"))
	}
	return(list(childTb=childTb[childIdx,, drop=F], parTb=parTb[parIdx,,drop=F]))
	
}


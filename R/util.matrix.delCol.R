util.matrix.delCol <-
function(ma, colPos){
	rnm = dim(ma)[1]
	cnm = dim(ma)[2]
	cutPt = (colPos-1)*rnm
	if(colPos==cnm){
	    re = matrix(ma[1:cutPt], nrow = rnm, ncol=(cnm-1))
	    return(re)
	}
        if(colPos == 1){
            re = matrix(ma[(rnm+1):(rnm*cnm)], nrow = rnm, ncol=(cnm-1))
	    return(re)
        }else{	
	    re = c(ma[1:cutPt], ma[(cutPt+rnm+1):(rnm*cnm)])
	    re = matrix(re, nrow = rnm, ncol=(cnm-1))									
	    return(re)
	}																			
}


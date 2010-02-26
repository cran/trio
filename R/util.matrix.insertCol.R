util.matrix.insertCol <-
function(ma, insertVec, afterCol){
	rnm = dim(ma)[1]
	cnm = dim(ma)[2]
	startnm = afterCol*rnm+1
        re = NULL
	if(afterCol==cnm){
	    re = matrix(c(ma, insertVec), ncol=cnm+1)
	}else{	
	    re = c(ma[1:startnm-1], insertVec, ma[startnm:(rnm*cnm)])
	    re = matrix(re, nrow = rnm, ncol=cnm+1)									
	}
        return(re)
      }


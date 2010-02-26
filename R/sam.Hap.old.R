sam.Hap.old <-
function(exp, prob, subjectCt){
	
	if(length(exp)==1){
		reExp = rep(exp, subjectCt)
	}else{
		reExp = sample(exp, size=subjectCt, prob=prob, replace=T)
	}
	
	return(reExp)
	
}


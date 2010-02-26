getP.logodds <-
function(logodds){
	return (exp(logodds)/(1+exp(logodds)))
}


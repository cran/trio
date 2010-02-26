sampleHapExclude <-
function(idxs, hapProb){
	newHapProb = hapProb
	newHapProb[idxs]=0
	idx = sample(1:length(newHapProb), prob=newHapProb, size=1)
	return(idx)
}


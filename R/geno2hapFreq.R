geno2hapFreq <-
function(prekey, block, genoFreq, snpBase=1, snpSeq){
	df = data.frame( key=I(rep(paste(prekey, block, sep="-"), 2)),
			hapF.ch = rep(prekey, 2),
			hapF.block = rep(block, 2),
			hapF.hap = c(snpBase, snpBase+1),
			freq = c(genoFreq[1]+genoFreq[2]/2, 1- (genoFreq[1]+genoFreq[2]/2) ),
			hapF.hapLen = rep(1, 2),
			hapF.markers_b = rep(snpSeq, 2),
			hapF.markers_e = rep(snpSeq, 2))
	
	return(df)
}


bindHapBkGenoMaps <-
function (hapBkMap=NULL, genoMap) 
{
  #print(str(genoMap))
  if(is.null(hapBkMap)){
		# when all the marker are single marker
		snpCtNum = nrow(genoMap$df)/3
		genomeMarkerInfo = matrix(NA, ncol=5, nrow=snpCtNum)
		genomeMarkerInfo[,1]=rep(1, snpCtNum)
		genomeMarkerInfo[,2]=1:snpCtNum
		genomeMarkerInfo[,3]=1:snpCtNum
		genomeMarkerInfo[,4]=1:snpCtNum
		genomeMarkerInfo[,5]=rep(1, snpCtNum)
                
	genoMapDf = cbind(genoMap$df, hapLens = rep(1, nrow(genoMap$df)))
	genoOnlyMap = dfToHapBkMap(genoMapDf, keyCol = NULL, chCol = genoMap$chCol, 
			blockCol = genoMap$genomeSeqCol, expCol = genoMap$expCol, 
			probCol = genoMap$probCol, hapLenCol = genoMap$genomeSeqCol, 
			beginCol = NULL, endCol = NULL, snpBase = genoMap$snpBase, 
			re.bf = T, re.javaGUI = T)
                
		genoIndex = which(genomeMarkerInfo[, 5] == 1)
		re = list(hapBkOnlyMap=NULL, genoOnlyMap = genoOnlyMap, genomeMarkerInfo = genomeMarkerInfo, genoIndex = genoIndex)
		return(re)
		
	}

	if (is.null(hapBkMap$beginCol)) {
		stop("hapBkMap$beginCol is null, methods not implemented.")
	}
	if (is.null(hapBkMap$dfStr)) {
		stop("hapBkMap$dfStr is null, methods not implemented.")
	}
	hapOnlyMap = hapBkMap
	genoMapDf = cbind(genoMap$df, hapLens = rep(1, nrow(genoMap$df)))
	genoOnlyMap = dfToHapBkMap(genoMapDf, keyCol = NULL, chCol = genoMap$chCol, 
			blockCol = genoMap$genomeSeqCol, expCol = genoMap$expCol, 
			probCol = genoMap$probCol, hapLenCol = genoMap$genomeSeqCol, 
			beginCol = NULL, endCol = NULL, snpBase = genoMap$snpBase, 
			re.bf = T, re.javaGUI = T)
	qu = NULL
	qu.type = NULL
	markers_gb = NULL
	markers_ge = NULL
	markersIndex_genome = 0
	oldCh = "ooo"
	newCh = "ooo"
	hapIdx = 0
	chkey = unlist(lapply(genoOnlyMap$keys, FUN = function(item) {
						re = util.str.seqCutter(item, delims = "-")[1]
						re
					}))
	chkey.hap = unlist(lapply(hapOnlyMap$keys, FUN = function(item) {
						re = util.str.seqCutter(item, delims = "-")[1]
						re
					}))
	chkeyUni = unique(chkey)
	chkeyUni.hap = unique(chkey.hap)
	chkeyDiff = setdiff(chkeyUni, chkeyUni.hap)
	if (length(chkeyDiff) == 0) {
		chkeyDiff = NULL
	}
	genomeMarkerInfo = NULL
	
	for (ikey in chkeyUni) {
		#print(ikey)
		snp.allCt = sum(chkey == ikey)
		if (is.element(ikey, chkeyDiff)) {
			tmp.endCt = sum(chkey == ikey)
			tmp.ma = NULL
			if (tmp.endCt > 0) {
				for (tt2 in 1:tmp.endCt) {
					tmp.ma = rbind(tmp.ma,
							c(rep(markersIndex_genome + tt2, 3), 1))
				}
				tmp.ma = matrix(tmp.ma, ncol = 4)
				hap.ma = data.frame(ch = I(rep(ikey, times = nrow(tmp.ma))), 
						markers_gb = tmp.ma[, 1], markers_ge = tmp.ma[, 1], 
						qu = tmp.ma[, 1], qu.type = tmp.ma[, 4])
			}
			else {
				stop("This is wrong.")
			}
			hap.ma = hap.ma[order(hap.ma[, 2]), ]
			genomeMarkerInfo = rbind(genomeMarkerInfo, hap.ma)
			markersIndex_genome = markersIndex_genome + snp.allCt
		}
		else {
			ibks.f <- (chkey.hap == ikey)
			
			tmp.ma = NULL
			## need to know whether there are some singleton appearing at the beginning
			if( hapOnlyMap$markers_b[ibks.f][1] > 1){
				for (tt in 1:(hapOnlyMap$markers_b[ibks.f][1]-1) ) {
					tmp.ma = rbind(tmp.ma, c(rep(markersIndex_genome+tt, 3), 1))
				}
			}
			
			markers_gb = markersIndex_genome + hapOnlyMap$markers_b[ibks.f]
			markers_ge = markersIndex_genome + hapOnlyMap$markers_e[ibks.f]
			qu = (hapIdx + 1):(hapIdx + sum(ibks.f))
			qu.type = rep(0, times = sum(ibks.f))
			hap.ma = data.frame(ch = I(rep(ikey, times = sum(ibks.f))), 
					markers_gb, markers_ge, qu, qu.type)
			hapIdx = hapIdx + sum(ibks.f)
			
			if (sum(ibks.f) > 1) {
				tmp.diff = markers_gb[-1] -
						markers_ge[-(sum(ibks.f))] - 1
				tmp.pos1 = which(tmp.diff > 0)
				if (length(tmp.pos1)>0) {
					tmp.pos = hap.ma[tmp.pos1, c(1, 3), drop = F]
					for (tt in 1:(length(tmp.pos1))) {
						for (tt2 in 1:(tmp.diff[tmp.pos1][tt])) {
							tmp.ma = rbind(tmp.ma, c(rep(tmp.pos[tt, 
															2] + tt2, 3), 1))
						}
					}
				}
			}
			tmp.endCt = snp.allCt - hapOnlyMap$markers_e[ibks.f][sum(ibks.f)]
			if (tmp.endCt > 0) {
				for (tt2 in 1:tmp.endCt) {
					tmp.ma = rbind(tmp.ma, c(
									rep(markers_ge[sum(ibks.f)] + 
													tt2, 3), 1))
				}
			}
			if (!is.null(tmp.ma)) {
				tmp.ma = matrix(tmp.ma, ncol = 4)
				single.ma = data.frame(ch = I(rep(ikey, times = nrow(tmp.ma))), 
						markers_gb = tmp.ma[, 1], markers_ge = tmp.ma[, 1], 
						qu = tmp.ma[, 1], qu.type = tmp.ma[, 4])
				hap.ma = rbind(hap.ma, single.ma)
				hap.ma = hap.ma[order(hap.ma[, 2]), ]
			}
			markersIndex_genome = markersIndex_genome + snp.allCt
			genomeMarkerInfo = rbind(genomeMarkerInfo, hap.ma)
		}
	}
	genomeMarkerInfo = genomeMarkerInfo[order(genomeMarkerInfo[, 2]), ]
	re = list(hapBkOnlyMap = hapOnlyMap, genoOnlyMap = genoOnlyMap, 
			genomeMarkerInfo = genomeMarkerInfo)
	hapIndex = which(genomeMarkerInfo[, 5] == 0)
	genoIndex = which(genomeMarkerInfo[, 5] == 1)
	re = c(re, list(hapIndex = hapIndex, genoIndex = genoIndex))
	return(re)
}


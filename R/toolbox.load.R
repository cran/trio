toolbox.load <-
function(freqMaps){

	## HARD CODE!!!HARD CODE assuming the genotype is never used in the genoMap
	genoMap = dfToGenoMap(df=freqMaps$genoMap.info, dataHeaders=c("prekey", "seq", "freq"),  genotype=c("11", "12", "22"), snpBase=1)
	
    ## need to build the overall map
    if(!is.null(freqMaps$hapMap.info)){
	  if(max(freqMaps$hapMap.info$hapLen>7)==1) stop("Cannot process haplotype block with 8 or more loci.")  
      bkFrame = freqMaps$hapMap.info
      bkFrame$key = paste(bkFrame[,1], bkFrame[,2], sep="-") 

      hapBkMap = dfToHapBkMap(data=bkFrame,  keyCol=ncol(bkFrame),
                   chCol=1, blockCol=2,  expCol=3, probCol=4, hapLenCol=5, beginCol=6, endCol=7, snpBase=1, re.bf = T, re.javaGUI = T)
      
	  hapBkGenoMap = bindHapBkGenoMaps(hapBkMap=hapBkMap, genoMap=genoMap)
		   
      hapBkMap = hapBkGenoMap$hapBkOnlyMap
	  changedKey =  hapBkMap$keys
		   
	  newHapBkMap = hapBkMap
	  for( i in changedKey){
			   newHapBkMap =  semiAugBkFrame(newHapBkMap, key=i, probLeftOver = .01)
	  }
	  semiAugHapBkGenoMap = hapBkGenoMap
	  semiAugHapBkGenoMap$hapBkOnlyMap = newHapBkMap
    }else{
	  hapBkGenoMap = bindHapBkGenoMaps(hapBkMap=NULL, genoMap=genoMap)
	  semiAugHapBkGenoMap = hapBkGenoMap
	}
   
    ## HARD CODE!!!HARD CODE
    maxSNP=7
    idx4hapDigitAll = exIdxFromHap(maxSNP)
    exhaustHapExpAll = exhaustHapExp(maxSNP, re.str=F)[[1]]
    
    maxRow = 4000
    maxTbl = matrix(NA, ncol =6, nrow = maxRow)
    maxMateTbl = matrix(NA, ncol =2, nrow = maxRow)
    
    varPrefix = paste(sample(letters[1:26], size=5), sep="", collapse="")
    varOriginalName = c("semiAugHapBkGenoMapSZ", "idx4hapDigitAll", "exhaustHapExpAll", "maxTbl", "maxMateTbl")
    appVarNames = paste(varPrefix, varOriginalName, sep="_")
    names(appVarNames) = c("freqMap", "digit", "exp", "tbl", "mateTbl")
    appVarNames = as.list(appVarNames)
    
    #print(paste("Application wise global environment var with names as", paste(appVarNames, sep="", collapse=";")))
    
    lapply(1:length(appVarNames), FUN=function(item, varNames, varList){
              assign(varNames[[item]], varList[[item]], env=.GlobalEnv); return(NULL)},
           varNames=appVarNames,
           varList = list(semiAugHapBkGenoMap, idx4hapDigitAll, exhaustHapExpAll, maxTbl, maxMateTbl))
    return(appVarNames)
}


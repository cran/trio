hapGenoBlockProc <-
function(hapGenoSub,  snpCoding=c(0, 1, 2, 3)){
	ifD = F
	if(ifD) print(hapGenoSub)
	
	snpCt = length(hapGenoSub)
	
	homoIn = NULL
	homoDigit = NULL
	missingIn = NULL
	hetoIn = NULL
	for ( i in 1: snpCt ){
		if(ifD) print(hapGenoSub[i])
		matchIndex = which(as.numeric(hapGenoSub[i])==snpCoding)
		if(matchIndex==1){
			missingIn = c(missingIn, i)
		}else if(matchIndex==2 | matchIndex==3){
			homoIn = c(homoIn, i)
			homoDigit = c(homoDigit, snpCoding[matchIndex])
		}else if(matchIndex==4){
			hetoIn = c(hetoIn, i)
		}  
	}
	return(list(homoIn=homoIn, homoDigit = homoDigit, missingIn=missingIn, hetoIn=hetoIn, ori=hapGenoSub))
}


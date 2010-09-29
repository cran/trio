# snp, snp1, snp2: vector of length 3*t, where t is the number of trios and each
#	of the n blocks of SNPs consists of the genotypes of father, mother, and
#	the offspring (in this order). Genotypes need to be coded by 0, 1, 2 (i.e.
#	the number of minor alleles). Missing values are allowed and need to be
#	coded by NA. Thus, the vector must have the same structure as the SNP columns
#	in the output of trio.check, or the genotype example data sets such as 
#	trio.gen1	


tdt <- function(snp, model=c("additive", "dominant", "recessive")){
	require(survival)
	n <- length(snp)
	if(n%%3 != 0)
		stop("snp has a length not dividable by 3.")
	if(any(!snp %in% c(0:2, NA)))
		stop("The values in snp must be 0, 1, and 2.")
	type <- match.arg(model)
	mat.trio <- matrix(snp, ncol=3, byrow=TRUE)
	mat.trio <- mat.trio[rowSums(is.na(mat.trio))==0,]
	mat <- matrix(c(0,0,0,0, 0,0,1,1, 0,0,1,1, 1,0,0,1, 1,0,0,1, 1,1,1,1, 1,1,1,1,
		2,2,2,2, 1,1,2,2, 1,1,2,2, 2,1,1,2, 2,1,1,2, 0,1,1,2, 1,0,1,2, 
		2,0,1,1), nrow=4)
	cn <- c("000", "010", "100", "011", "101", "021", "201", "222", "121", "211",
		"122", "212", "110", "111", "112")
	colnames(mat) <- cn
	code <- paste(mat.trio[,1], mat.trio[,2], mat.trio[,3], sep="")
	if(any(!code %in% cn)){
		tmp.ids <- !code %in% cn
		warning(sum(tmp.ids), " trios show Mendelian errors. These are removed.",
			call.=FALSE)
		code <- code[!tmp.ids]
	}
	x <- as.vector(mat[,code])
	if(type=="dominant")
		x <- (x >= 1) * 1
	if(type=="recessive")
		x <- (x > 1) * 1
	y <- rep.int(c(1,0,0,0), length(code))
	strat <- rep(1:length(code), e=4)
	c.out <- clogit(y ~ x + strata(strat))
	coef <- c.out$coefficients
	names(coef) <- NULL
	se <- sqrt(diag(c.out$var))
	stat <- (coef/se)^2
	conf <- exp(coef + c(-1, 1) * qnorm(0.975) * se)
	out <- list(coef=coef, se=se, stat=stat, pval= 1 - pchisq(stat, 1),
		OR=exp(coef), lowerOR=conf[1], upperOR=conf[2], ia=FALSE, type=type)  
	class(out) <- "tdt"
	out
}

print.tdt <- function(x, digits=4, ...){
	pval <- format.pval(x$pval, digits=digits)
	out <- data.frame(Coef=x$coef, OR=x$OR, Lower=x$lowerOR, Upper=x$upperOR, 
		SE=x$se, Statistic=x$stat, "p-Value"=pval, check.names=FALSE)
	if(length(x$coef)==1)
		rownames(out) <- ""
	if(x$ia)
		cat("Genotypic TDT for Two-Way Interaction (Using 15 Pseudo Controls)",
			"\n\n")
	else
		cat("        Genotypic TDT Based on 3 Pseudo Controls","\n\n")
	cat("Model Type:", switch(x$type, "additive"="Additive", "dominant"="Dominant",
		"recessive"="Recessive"), "\n\n")
	if(x$ia && is.na(x$stat))
		cat("Fitting failed. Possible reason: SNPs might be in (strong) LD.\n\n")
	
	else
		print(format(out, digits=digits))
	if(!is.na(x$pval) && x$pval <= .Machine$double.eps)
		warning("The results might be misleading, as the very small p-value might be caused\n",
			"by non-biological artefacts (e.g., sparseness of the data).", call.=FALSE)
}

tdt2way <- function(snp1, snp2, epistatic=TRUE,
		model=c("additive", "dominant", "recessive"), add=FALSE,
		warnError=TRUE){
	require(survival)
	n1 <- length(snp1)
	n2 <- length(snp2)
	if(n1 != n2)
		stop("snp1 must have the same length as snp2.")
	if(n1%%3 != 0)
		stop("The SNPs do not seem to contain trio data, as their length\n",
			"is not dividable by 3.")
	if(any(!snp1 %in% c(0, 1, 2, NA)))
		stop("The values in snp1 must be 0, 1, and 2.")
	if(any(!snp2 %in% c(0, 1, 2, NA)))
		stop("The values in snp2 must be 0, 1, and 2.")
	type <- match.arg(model)
	mat.trio <- cbind(matrix(snp1, ncol=3, byrow=TRUE), matrix(snp2, ncol=3, byrow=TRUE))
	mat.trio <- mat.trio[rowSums(is.na(mat.trio))==0,]
	mat <- matrix(c(0,0,0,0, 0,0,1,1, 0,0,1,1, 1,0,0,1, 1,0,0,1, 1,1,1,1, 1,1,1,1,
		2,2,2,2, 1,1,2,2, 1,1,2,2, 2,1,1,2, 2,1,1,2, 0,1,1,2, 1,0,1,2, 
		2,0,1,1), nrow=4)
	cn <- c("000", "010", "100", "011", "101", "021", "201", "222", "121", "211",
		"122", "212", "110", "111", "112")
	colnames(mat) <- cn
	code1 <- paste(mat.trio[,1], mat.trio[,2], mat.trio[,3], sep="")
	code2 <- paste(mat.trio[,4], mat.trio[,5], mat.trio[,6], sep="")
	if(any(!code1 %in% cn)){
		tmp.ids <- !code1 %in% cn
		warning(sum(tmp.ids), " trios in snp1 show Mendelian errors. These are removed.",
			call.=FALSE)
		code1 <- code1[!tmp.ids]
		code2 <- code2[!tmp.ids]
	}
	if(any(!code2 %in% cn)){
		tmp.ids <- !code2 %in% cn
		warning(sum(tmp.ids), " trios show Mendelian errors in snp2 (but not in snp1).\n",
			"These trios are removed.", call.=FALSE)
		code2 <- code2[!tmp.ids]
		code1 <- code1[!tmp.ids]
	}  
	n.trio <- length(code1)		
	x1 <- as.vector(mat[,code1])
	x2 <- as.vector(mat[,code2])
	if(epistatic)
		return(tdtEpistatic(x1, x2, n.trio, warnError=warnError))
	if(type=="dominant"){
		x1 <- (x1 >= 1) * 1
		x2 <- (x2 >= 1) * 1
	}
	if(type=="recessive"){
		x1 <- (x1 > 1) * 1
		x2 <- (x2 > 1) * 1
	} 
	mat.ids <- cbind(rep.int(rep(1:4, e=4), n.trio), rep.int(1:4, 4*n.trio))
	mat.ids <- mat.ids + 4 * rep(0:(n.trio-1), e=16)
	IA <- x1[mat.ids[,1]] * x2[mat.ids[,2]]
	y <- rep.int(c(1, rep.int(0, 15)), n.trio)
	strat <- rep(1:n.trio, e=16)
	if(warnError){
		wa <- options()$warn
		on.exit(options(warn=wa))
		options(warn=2)
	}
	if(!add)
		c.out <- try(clogit(y ~ IA + strata(strat)), silent=TRUE)
	else{
		SNP1 <- x1[mat.ids[,1]]
		SNP2 <- x2[mat.ids[,2]]
		c.out <- try(clogit(y ~ SNP1 + SNP2 + IA + strata(strat)), silent=TRUE)
	}
	if(warnError)
		options(warn=wa)
	if(is(c.out, "try-error"))
		coef <- se <- stat <- pval <- lower <- upper <- NA
	else{
		coef <- c.out$coefficients
		se <- sqrt(diag(c.out$var))
		if(!add)
			names(coef) <- NULL
		stat <- (coef/se)^2
		pval <- pchisq(stat, 1, lower.tail=FALSE)
		lower <- exp(coef - qnorm(0.975) * se)
		upper <- exp(coef + qnorm(0.975) * se)
	}
	out <- list(coef=coef, se=se, stat=stat, pval= pval,
		OR=exp(coef), lowerOR=lower, upperOR=upper, ia=TRUE, type=type)  
	class(out) <- "tdt"
	out

}
	

colTDT <- function(mat.snp, model=c("additive", "dominant", "recessive"), fast=TRUE, size=50,  
		warnError=TRUE){
	if(!is.matrix(mat.snp))
		stop("mat.snp has to be a matrix.")
	if(nrow(mat.snp) %% 3 != 0)
		stop("mat.snp does not seem to contain trio data, as its number of rows is\n",
			"not dividable by 3.")
	if(is.null(rownames(mat.snp)))
		stop("mat.snp does not seem to be a matrix in genotype format,\n",
			"as the names of the rows are missing.")
	if(any(!mat.snp %in% c(0, 1, 2, NA)))
		stop("The values in mat.snp must be 0, 1, and 2.")
	type <- match.arg(model)
	if(fast)
		return(fastTDT(mat.snp, type, size=size))	
	require(survival)
	warning("This way (fast = FALSE) of performing a gTDT will not be available in\n",
		"future versions of the package trio.")
	mat.code <- matrix(c(0,0,0,0, 0,0,1,1, 0,0,1,1, 1,0,0,1, 1,0,0,1, 1,1,1,1, 
		1,1,1,1, 2,2,2,2, 1,1,2,2, 1,1,2,2, 2,1,1,2, 2,1,1,2, 0,1,1,2, 
		1,0,1,2, 2,0,1,1), nrow=4)
	cn <- c("000", "010", "100", "011", "101", "021", "201", "222", "121", "211",
		"122", "212", "110", "111", "112")
	colnames(mat.code) <- cn
	coef <- se <- numeric(ncol(mat.snp))
	if(warnError){
		wa <- options()$warn
		on.exit(options(warn=wa))
	}
	for(i in 1:ncol(mat.snp)){
		mat.trio <- matrix(mat.snp[,i], ncol=3, byrow=TRUE)
		mat.trio <- mat.trio[rowSums(is.na(mat.trio))==0,]
		code <- paste(mat.trio[,1], mat.trio[,2], mat.trio[,3], sep="")
		if(any(!code %in% cn)){
			tmp.ids <- !code %in% cn
			warning(sum(tmp.ids), " trios show Mendelian errors. These are removed.",
				call.=FALSE)
			code <- code[!tmp.ids]
		}
		x <- as.vector(mat.code[,code])
		if(type=="dominant")
			x <- (x >= 1) * 1
		if(type=="recessive")
			x <- (x > 1) * 1
		y <- rep.int(c(1,0,0,0), length(code))
		strat <- rep(1:length(code), e=4)
		if(warnError)
			options(warn=2)
		c.out <- try(clogit(y ~ x + strata(strat)), silent=TRUE)
		if(is(c.out, "try-error"))
			coef[i] <- se[i] <- NA
		else{
			coef[i] <- c.out$coefficients
			se[i] <- sqrt(diag(c.out$var))
		}
		if(warnError)
			options(warn=wa)
	}
	if(any(is.na(coef)))
		warning("The fitting of some of the models has failed. A possible reason\n",
			"is that the corresponding SNPs have very low minor allele frequencies.\n",
			"For these SNPs, all statistics are thus set to NA.", call.=FALSE)
	stat <- (coef / se)^2
	lower <- exp(coef - qnorm(0.975) * se)
	upper <- exp(coef + qnorm(0.975) * se)
	pval <- 1 - pchisq(stat, 1)
	if(is.null(colnames(mat.snp)))
		names(coef) <- names(pval) <- names(stat) <- paste("SNP", 1:ncol(mat.snp), sep="")
	else
		names(coef) <- names(pval) <- names(stat) <- colnames(mat.snp)
	out <- list(coef=coef, se=se, stat=stat, pval=pval, OR=exp(coef), 
		lowerOR=lower, upperOR=upper, ia=FALSE, type=type, add=FALSE) 
	class(out) <- "colTDT"
	out
}

colTDT2way <- function(mat.snp, epistatic=TRUE, genes=NULL, maf=FALSE,
		model=c("additive", "dominant", "recessive"), add=FALSE,
		warnError=TRUE){
	require(survival)
	if(!is.matrix(mat.snp))
		stop("mat.snp has to be a matrix.")
	if(nrow(mat.snp) %% 3 != 0)
		stop("mat.snp does not seem to contain trio data, as its number of rows is\n",
			"not dividable by 3.")
	if(is.null(rownames(mat.snp)))
		stop("mat.snp does not seem to be a matrix in genotype format,\n",
			"as the names of the rows are missing.")
	if(any(!mat.snp %in% c(0, 1, 2, NA)))
		stop("The values in mat.snp must be 0, 1, and 2.")
	type <- match.arg(model)
	mat.code <- matrix(c(0,0,0,0, 0,0,1,1, 0,0,1,1, 1,0,0,1, 1,0,0,1, 1,1,1,1, 
		1,1,1,1, 2,2,2,2, 1,1,2,2, 1,1,2,2, 2,1,1,2, 2,1,1,2, 0,1,1,2, 
		1,0,1,2, 2,0,1,1, NA,NA,NA,NA), nrow=4)
	cn <- c("000", "010", "100", "011", "101", "021", "201", "222", "121", "211",
		"122", "212", "110", "111", "112", "NANANA")
	colnames(mat.code) <- cn
	n.snp <- ncol(mat.snp)
	mat.pseudo <- matrix(NA, 4/3 * nrow(mat.snp), n.snp)
	for(i in 1:n.snp){
		mat.trio <- matrix(mat.snp[,i], ncol=3, byrow=TRUE)
		mat.trio[rowSums(is.na(mat.trio)) > 0, ] <- NA
		code <- paste(mat.trio[,1], mat.trio[,2], mat.trio[,3], sep="")
		if(any(!code %in% cn)){
			tmp.ids <- !code %in% cn
			warning(sum(tmp.ids), " trios show Mendelian errors. These are removed.",
				call.=FALSE)
			code[tmp.ids] <- "NANANA"
		}
		mat.pseudo[,i] <- as.vector(mat.code[,code])
	}
	if(maf){
		mat.snp <- mat.snp[-seq(3, nrow(mat.snp), 3), ]
		mat.snp <- mat.snp[!duplicated(rownames(mat.snp)), ]
		valMAF <- colSums(mat.snp, na.rm=TRUE) / (2 * colSums(!is.na(mat.snp)))
	}
	else
		valMAF <- NULL
	n.trio <- length(code)
	if(epistatic)
		return(colTDTepistatic(mat.pseudo, n.trio, colnames(mat.snp), genes=genes, valMAF=valMAF,
			warnError=warnError))
	if(type=="dominant")
		mat.pseudo <- (mat.pseudo >= 1) * 1
	if(type=="recessive")
		mat.pseudo <- (mat.pseudo > 1) * 1
	mat.ids <- cbind(rep.int(rep(1:4, e=4), n.trio), rep.int(1:4, 4*n.trio))
	mat.ids <- mat.ids + 4 * rep(0:(n.trio-1), e=16)
	y <- rep.int(c(1, rep.int(0, 15)), n.trio)
	strat <- rep(1:n.trio, e=16)
	if(is.null(genes)){
		coef <- se <- numeric(n.snp * (n.snp - 1) / 2)
		combs <- allCombs(n.snp)
	}
	else{
		if(!is.character(genes))
			stop("genes must be a vector of character strings.")
		if(length(genes) != n.snp)
			stop("The length of genes must be equal to the number of columns of mat.snp.")
		ids.genes <- as.numeric(as.factor(genes))
		combs <- allBetweenCombs(ids.genes)
		coef <- se <- numeric(nrow(combs))
	}
	if(warnError){
		wa <- options()$warn
		on.exit(options(warn=wa))
		options(warn=2)
	}
	for(i in 1:nrow(combs)){
		x1 <- mat.pseudo[mat.ids[,1], combs[i,1]] 
		x2 <- mat.pseudo[mat.ids[,2], combs[i,2]]
		IA <- x1 * x2
		if(add)
			c.out <- try(clogit(y ~ IA + x1 + x2 + strata(strat)), silent=TRUE)
		else
			c.out <- try(clogit(y ~ IA + strata(strat)), silent=TRUE)
		coef[i] <- if(is(c.out, "try-error")) NA else c.out$coefficients[1]
		se[i] <- if(is(c.out, "try-error")) NA else sqrt(diag(c.out$var))[1]
	}
	if(warnError)
		options(warn=wa)
	if(any(is.na(coef)))
		warning("The fitting of some of the models has failed. A possible reason\n",
			"is that the two SNPs might be in (strong) LD. For these interactions,\n",
			"all statistics are therefore set to NA.", call.=FALSE)
	stat <- (coef / se)^2
	lower <- exp(coef - qnorm(0.975) * se)
	upper <- exp(coef + qnorm(0.975) * se)
	pval <- pchisq(stat, 1, lower.tail=FALSE)
	if(any(pval <= .Machine$double.eps, na.rm=TRUE))
		warning("Some of the p-values are very small (< 2.2e-16). These results might be\n",
			"misleading, as the small p-values might be caused by non-biological artefacts\n",
			"(e.g., sparseness of the data).", call.=FALSE)
	rn <- if(is.null(colnames(mat.snp))) paste("SNP", 1:n.snp, sep="") else colnames(mat.snp)
	names(coef) <- names(stat) <- names(pval) <- paste(rn[combs[,1]], rn[combs[,2]], sep=" : ")
	if(!is.null(genes))
		genes <- paste(genes[combs[,1]], genes[combs[,2]], sep=" : ")
	if(maf){
		mat.maf <- cbind(round(valMAF, 4)[combs[,1]], round(valMAF, 4)[combs[,2]])
		rownames(mat.maf) <- names(coef) 
		colnames(mat.maf) <- c("First SNP", "Second SNP")
	}
	else
		mat.maf <- NULL
	out <- list(coef=coef, se=se, stat=stat, pval=pval, OR=exp(coef), 
		lowerOR=lower, upperOR=upper, ia=TRUE, type=type, add=add, genes=genes, maf=valMAF,
		matMAF=mat.maf) 
	class(out) <- "colTDT"
	out	

}

print.colTDT <- function(x, top=5, digits=4, ...){
	pval <- format.pval(x$pval, digits=digits)
	out <- data.frame(Coef=x$coef, OR=x$OR, Lower=x$lowerOR, Upper=x$upperOR, 
		SE=x$se, Statistic=x$stat, "p-Value"=pval, check.names=FALSE)
	if(x$ia)
		cat("      Genotypic TDT for Two-Way Interaction (Using 15 Pseudo Controls)",
			"\n\n")
	else
		cat("        Genotypic TDT Based on 3 Pseudo Controls","\n\n")
	cat("Model Type: ", switch(x$type, "additive"="Additive", "dominant"="Dominant",
		"recessive"="Recessive"), "\n", 
		if(x$add) "Model also contains the two respective individual SNPs.\n",
		"\n", sep="")
	if(length(x$coef) > top){
		ord <- order(x$pval)[1:top]
		out <- out[ord,]
		cat("Top", top, ifelse(x$ia, "SNP Interactions:\n", "SNPs:\n"))
	}
	print(format(out, digits=digits))
}

tdtEpistatic <- function(x1, x2, n.trio, warnError=TRUE){
	x1 <- x1 - 1
	z1 <- ifelse(x1==0, 0.5, -0.5)
	x2 <- x2 - 1
	z2 <- ifelse(x2==0, 0.5, -0.5)

	mat.ids <- cbind(rep.int(rep(1:4, e=4), n.trio), rep.int(1:4, 4*n.trio))
	mat.ids <- mat.ids + 4 * rep(0:(n.trio-1), e=16)
	x1 <- x1[mat.ids[,1]]
	x2 <- x2[mat.ids[,2]]
	z1 <- z1[mat.ids[,1]]
	z2 <- z2[mat.ids[,2]]
	x1x2 <- x1 * x2
	x1z2 <- x1 * z2
	z1x2 <- z1 * x2
	z1z2 <- z1 * z2
	y <- rep.int(c(1, rep.int(0, 15)), n.trio)
	strat <- rep(1:n.trio, e=16)
	wa <- options()$warn
	on.exit(options(warn=wa))
	options(warn=ifelse(warnError, 2, wa))
	woIA <- try(clogit(y ~ x1 + z1 + x2 + z2 + strata(strat)), silent=TRUE)
	full <- try(clogit(y ~ x1 + z1 + x2 + z2 + x1x2 + x1z2 + z1x2 + z1z2 + 
		strata(strat)), silent=TRUE)
	options(warn=wa)
	ll.main <- if(is(woIA, "try-error")) NA else woIA$loglik[2]	
	ll.full <- if(is(full, "try-error")) NA else full$loglik[2]
	if(!is.na(ll.full) && !is.na(ll.main)){
		stat <- -2 * (ll.main - ll.full)
		pval <- pchisq(stat, 4, lower.tail=FALSE)
	}
	else
		stat <- pval <- NA
	out <- list(ll.main=ll.main, ll.full=ll.full, stat=stat, pval=pval,
		full.model=full) 
	class(out) <- "tdtEpi"
	out
}

colTDTepistatic <- function(mat.pseudo, n.trio, rn, genes=NULL, valMAF=NULL, warnError=TRUE){
	x <- mat.pseudo - 1 
	z <- (x == 0) - 0.5
	mat.ids <- cbind(rep.int(rep(1:4, e=4), n.trio), rep.int(1:4, 4*n.trio))
	mat.ids <- mat.ids + 4 * rep(0:(n.trio-1), e=16)
	y <- rep.int(c(1, rep.int(0, 15)), n.trio)
	strat <- rep(1:n.trio, e=16)
	n.snp <- ncol(mat.pseudo)
	if(is.null(genes)){
		ll.main <- ll.full <- numeric(n.snp * (n.snp - 1) / 2)
		combs <- allCombs(n.snp)
	}
	else{
		if(!is.character(genes))
			stop("genes must be a vector of character strings.")
		if(length(genes) != n.snp)
			stop("The length of genes must be equal to the number of columns of mat.snp.")
		ids.genes <- as.numeric(as.factor(genes))
		combs <- allBetweenCombs(ids.genes)
		ll.main <- ll.full <- numeric(nrow(combs))
	}
	if(warnError){
		wa <- options()$warn
		on.exit(options(warn=wa))
		options(warn=2)
	}
	for(i in 1:nrow(combs)){
		x1 <- x[mat.ids[,1], combs[i,1]]
		z1 <- z[mat.ids[,1], combs[i,1]] 
		x2 <- x[mat.ids[,2], combs[i,2]]
		z2 <- z[mat.ids[,2], combs[i,2]]
		woIA <- try(clogit(y ~ x1 + z1 + x2 + z2 + strata(strat)), silent=TRUE)
		full <- try(clogit(y ~ x1 + z1 + x2 + z2 + x1*x2 + x1*z2 + z1*x2 + 
			z1*z2 + strata(strat)), silent=TRUE)
		ll.main[i] <- if(is(woIA, "try-error")) NA else woIA$loglik[2]
		ll.full[i] <- if(is(full, "try-error")) NA else full$loglik[2]
	}
	if(warnError)
		options(warn=wa)
	if(any(is.na(ll.full)) | any(is.na(ll.main)))
		warning("The fitting of some of the models has failed. A possible reason\n",
			"is that the two respective SNPs might be in (strong) LD.\n",
			"The corresponding test statistic and the p-value are thus set to NA.",
			call.=FALSE)
	stat <- -2 * (ll.main - ll.full)
	pval <- pchisq(stat, 4, lower.tail=FALSE)
	if(any(pval <= .Machine$double.eps, na.rm=TRUE))
		warning("Some of the p-values are very small (< 2.2e-16). These results might be\n",
			"misleading, as the small p-values might be caused by non-biological artefacts\n",
			"(e.g., sparseness of the data).", call.=FALSE)
	if(is.null(rn))
		rn <- paste("SNP", 1:n.snp, sep="")
	names(ll.full) <- names(stat) <- names(pval) <- paste(rn[combs[,1]], rn[combs[,2]], sep=" : ")
	if(!is.null(genes)){
		ind.genes <- genes
		genes <- paste(genes[combs[,1]], genes[combs[,2]], sep=" : ")
	}
	else
		ind.genes <- NULL
	if(!is.null(valMAF)){
		mat.maf <- cbind(round(valMAF, 4)[combs[,1]], round(valMAF, 4)[combs[,2]])
		rownames(mat.maf) <- names(ll.full)
		colnames(mat.maf) <- c("First MAF", "Second MAF") 
	}
	else
		mat.maf <- NULL
	out <- list(ll.main=ll.main, ll.full=ll.full, stat=stat, pval=pval, 
		maf=valMAF, matMAF=mat.maf, genes=genes, ind.genes=ind.genes)
	class(out) <- "colTDTepi"
	out
}
	

print.tdtEpi <- function(x, digits=3, ...){
	pval <- format.pval(x$pval, digits=digits-1)
	cat("Genotypic TDT for Epistatic Interactions (Using 15 Pseudo Controls)", 
		"\n\n", sep="")
	cat("Likelihood Ratio Test:\n", sep="")
	if(is.na(x$stat))
		cat("Failed. A possible reason is that the SNPs might be in (strong) LD.",
			"\n\n")
	else{
		cat("Loglikelihood (with Interactions): ", round(x$ll.full, digits), "\n", sep="")
		cat("Loglikelihood (without IAs):       ", round(x$ll.main, digits), "\n", sep="")
		cat("Test Statistic: ", round(x$stat, digits), "\n", sep="")
		cat("P-Value:        ", pval, "\n\n", sep="")
	}
}

print.colTDTepi <- function(x, top=5, digits=4, ...){
	pval <- format.pval(x$pval, digits=digits)
	out <- data.frame("LL (with IAs)"=x$ll.full, "LL (w/o IAs)"=x$ll.main,
		Statistic=x$stat, "P-Value"=pval, check.names=FALSE)
	if(!is.null(x$genes))
		out <- data.frame(out, Genes=x$genes, check.names=FALSE, stringsAsFactors=FALSE)
	if(!is.null(x$matMAF))
		out <- data.frame(out, x$matMAF, check.names=FALSE, stringsAsFactors=FALSE) 
	cat("      Genotypic TDT for Epistatic Interactions (Using 15 Pseudo Controls)",
		"\n\n")
	if(length(x$ll.main) > top){
		ord <- order(x$pval)[1:top]
		out <- out[ord,]
		cat("Top", top, "SNP Interactions (Likelihood Ratio Test):\n")
	}	
	else
		cat("Likelihood Ratio Test:\n")
	print(format(out, digits=digits))
}
	


colTDTinter2way <- function(mat.snp1, mat.snp2, snp=NULL, epistatic=TRUE, maf=FALSE,
		model=c("additive", "dominant", "recessive"), add=FALSE,
		warnError=TRUE){
	require(survival)
	if(missing(mat.snp2)){
		if(is.null(snp))
			stop("Either mat.snp2 or snp has to be specified.")
		mat.snp2 <- cbind(snp)
	}
	else{
		if(!is.null(snp))
			stop("If mat.snp2 is specified, snp cannot be specified.")
	}
	if(!is.matrix(mat.snp1))
		stop("mat.snp1 has to be a matrix.")
	if(!is.matrix(mat.snp2))
		stop("mat.snp2 has to be a matrix.")
	n.obs <- nrow(mat.snp1)
	if(n.obs != nrow(mat.snp2))
		stop("mat.snp1 must have the same number of rows as mat.snp2.")
	if(n.obs %% 3 != 0)
		stop("The matrices do not seem to contain trio data, as its number of rows is\n",
			"not dividable by 3.")
	if(is.null(rownames(mat.snp1)))
		stop("mat.snp1 does not seem to be a matrix in genotype format,\n",
			"as the names of the rows are missing.")
	if(!is.null(rownames(mat.snp1)) && !is.null(rownames(mat.snp2)) && 
		any(rownames(mat.snp1) != rownames(mat.snp2)))
		stop("mat.snp1 must have the same row names as mat.snp2.")
	if(any(!mat.snp1 %in% c(0, 1, 2, NA)))
		stop("The values in mat.snp1 must be 0, 1, and 2.")
	if(any(!mat.snp2 %in% c(0, 1, 2, NA)))
		stop("The values in mat.snp2 must be 0, 1, and 2.")
	type <- match.arg(model)
	mat.code <- matrix(c(0,0,0,0, 0,0,1,1, 0,0,1,1, 1,0,0,1, 1,0,0,1, 1,1,1,1, 
		1,1,1,1, 2,2,2,2, 1,1,2,2, 1,1,2,2, 2,1,1,2, 2,1,1,2, 0,1,1,2, 
		1,0,1,2, 2,0,1,1, NA,NA,NA,NA), nrow=4)
	cn <- c("000", "010", "100", "011", "101", "021", "201", "222", "121", "211",
		"122", "212", "110", "111", "112", "NANANA")
	colnames(mat.code) <- cn
	n.snp1 <- ncol(mat.snp1)
	mat.pseudo1 <- matrix(NA, 4/3 * n.obs, n.snp1)
	for(i in 1:n.snp1){
		mat.trio <- matrix(mat.snp1[,i], ncol=3, byrow=TRUE)
		mat.trio[rowSums(is.na(mat.trio)) > 0, ] <- NA
		code <- paste(mat.trio[,1], mat.trio[,2], mat.trio[,3], sep="")
		if(any(!code %in% cn)){
			tmp.ids <- !code %in% cn
			warning(sum(tmp.ids), " trios in mat.snp1 show Mendelian errors. These are removed.",
				call.=FALSE)
			code[tmp.ids] <- "NANANA"
		}
		mat.pseudo1[,i] <- as.vector(mat.code[,code])
	}
	n.snp2 <- ncol(mat.snp2)
	mat.pseudo2 <- matrix(NA, 4/3 * n.obs, n.snp2)
	for(i in 1:n.snp2){
		mat.trio <- matrix(mat.snp2[,i], ncol=3, byrow=TRUE)
		mat.trio[rowSums(is.na(mat.trio)) > 0, ] <- NA
		code <- paste(mat.trio[,1], mat.trio[,2], mat.trio[,3], sep="")
		if(any(!code %in% cn)){
			tmp.ids <- !code %in% cn
			warning(sum(tmp.ids), " trios show Mendelian errors. These are removed.",
				call.=FALSE)
			code[tmp.ids] <- "NANANA"
		}
		mat.pseudo2[,i] <- as.vector(mat.code[,code])
	}
	if(maf){
		mat.tmp <- mat.snp1[-seq(3, nrow(mat.snp1), 3), ]
		mat.tmp <- mat.tmp[!duplicated(rownames(mat.tmp)),]
		valMAF1 <- colSums(mat.tmp, na.rm=TRUE) / (2 * colSums(!is.na(mat.tmp)))
		mat.tmp <- mat.snp2[-seq(3, nrow(mat.snp2), 3), ]
		mat.tmp <- mat.tmp[!duplicated(rownames(mat.tmp)),]
		valMAF2 <- colSums(mat.tmp, na.rm=TRUE) / (2 * colSums(!is.na(mat.tmp)))
		valMAF <- list(mat.snp1=valMAF1, mat.snp2=valMAF2)
	}
	else
		valMAF <- NULL
	n.trio <- length(code)
	if(epistatic)
		return(colTDTinterEpi(mat.pseudo1, mat.pseudo2, n.trio, colnames(mat.snp1), colnames(mat.snp2),
			valMAF=valMAF, warnError=warnError))
	if(type=="dominant"){
		mat.pseudo1 <- (mat.pseudo1 >= 1) * 1
		mat.pseudo2 <- (mat.pseudo2 >= 1) * 1
	}
	if(type=="recessive"){
		mat.pseudo1 <- (mat.pseudo1 > 1) * 1
		mat.pseudo2 <- (mat.pseudo2 > 1) * 1
	}
	mat.ids <- cbind(rep.int(rep(1:4, e=4), n.trio), rep.int(1:4, 4*n.trio))
	mat.ids <- mat.ids + 4 * rep(0:(n.trio-1), e=16)
	y <- rep.int(c(1, rep.int(0, 15)), n.trio)
	strat <- rep(1:n.trio, e=16)
	coef <- se <- numeric(n.snp1 * n.snp2)
	if(warnError){
		wa <- options()$warn
		on.exit(options(warn=wa))
		options(warn=2)
	}
	for(i in 1:n.snp1){
		x1 <- mat.pseudo1[mat.ids[,1], i]
		for(j in 1:n.snp2){ 
			x2 <- mat.pseudo2[mat.ids[,2], j]
			IA <- x1 * x2
			if(add)
				c.out <- try(clogit(y ~ IA + x1 + x2 + strata(strat)), silent=TRUE)
			else
				c.out <- try(clogit(y ~ IA + strata(strat)), silent=TRUE)
			k <- (i-1) * n.snp2 + j
			coef[k] <- if(is(c.out, "try-error")) NA else c.out$coefficients[1]
			se[k] <- if(is(c.out, "try-error")) NA else sqrt(diag(c.out$var))[1]
		}
	}
	if(warnError)
		options(warn=wa)
	if(any(is.na(coef)))
		warning("The fitting of some of the models has failed. A possible reason\n",
			"is that the two corresponding SNPs might be in (strong) LD.\n",
			"For these interactions, all statistics are therefore set to NA.", 
			call.=FALSE)
	stat <- (coef / se)^2
	lower <- exp(coef - qnorm(0.975) * se)
	upper <- exp(coef + qnorm(0.975) * se)
	pval <- pchisq(stat, 1, lower.tail=FALSE)
	if(any(pval <= .Machine$double.eps, na.rm=TRUE))
		warning("Some of the p-values are very small (< 2.2e-16). These results might be\n",
			"misleading, as the small p-values might be caused by non-biological artefacts\n",
			"(e.g., sparseness of the data).", call.=FALSE)
	rn1 <- if(is.null(colnames(mat.snp1))) paste("SNP", 1:n.snp1, sep="") else colnames(mat.snp1)
	if(is.null(snp)){
		rn2 <- if(is.null(colnames(mat.snp2))) paste("SNP", n.snp1 + (1:n.snp2), sep="") else colnames(mat.snp2)
		names(coef) <- names(stat) <- names(pval) <- paste(rep(rn1, e=n.snp2), rep(rn2, n.snp1), sep=" : ")
	}
	else
		names(coef) <- names(stat) <- names(pval) <- rn1
	if(maf){
		mat.maf <- cbind(rep(round(valMAF1, 4), e=n.snp2), rep(round(valMAF2, 4), n.snp1))
		rownames(mat.maf) <- names(coef)
		colnames(mat.maf) <- c("First SNP", "Second SNP")
	}
	else
		mat.maf <- NULL
	out <- list(coef=coef, se=se, stat=stat, pval=pval, OR=exp(coef), 
		lowerOR=lower, upperOR=upper, ia=TRUE, type=type, add=add,
		maf=valMAF, matMAF=mat.maf) 
	class(out) <- "colTDT"
	out	

}
	
colTDTinterEpi <- function(mat.pseudo1, mat.pseudo2, n.trio, rn1, rn2, valMAF=NULL,
		warnError=TRUE){
	mat.x1 <- mat.pseudo1 - 1 
	mat.z1 <- (mat.x1 == 0) - 0.5
	mat.x2 <- mat.pseudo2 - 1
	mat.z2 <- (mat.x2 == 0) - 0.5
	mat.ids <- cbind(rep.int(rep(1:4, e=4), n.trio), rep.int(1:4, 4*n.trio))
	mat.ids <- mat.ids + 4 * rep(0:(n.trio-1), e=16)
	y <- rep.int(c(1, rep.int(0, 15)), n.trio)
	strat <- rep(1:n.trio, e=16)
	n.snp1 <- ncol(mat.pseudo1)
	n.snp2 <- ncol(mat.pseudo2)
	ll.main <- ll.full <- numeric(n.snp1 * n.snp2)
	if(warnError){
		wa <- options()$warn
		on.exit(options(warn=wa))
		options(warn=2)
	}
	for(i in 1:n.snp1){
		x1 <- mat.x1[mat.ids[,1], i]
		z1 <- mat.z1[mat.ids[,1], i]
		for(j in 1:n.snp2){ 
			x2 <- mat.x2[mat.ids[,2], j]
			z2 <- mat.z2[mat.ids[,2], j]
			woIA <- try(clogit(y ~ x1 + z1 + x2 + z2 + strata(strat)), silent=TRUE)
			full <- try(clogit(y ~ x1 + z1 + x2 + z2 + x1*x2 + x1*z2 + z1*x2 + 
				z1*z2 + strata(strat)), silent=TRUE)
			k <- (i-1) * n.snp2 + j
			ll.main[k] <- if(is(woIA, "try-error")) NA else woIA$loglik[2]
			ll.full[k] <- if(is(full, "try-error")) NA else full$loglik[2]
		}
	}
	if(warnError)
		options(warn=wa)
	if(any(is.na(ll.full)) | any(is.na(ll.main)))
		warning("The fitting of some of the models has failed. A possible reason\n",
			"is that the two respective SNPs might be in (strong) LD.\n",
			"The corresponding test statistic and the p-value are thus set to NA.",
			call.=FALSE)
	stat <- -2 * (ll.main - ll.full)
	pval <- pchisq(stat, 4, lower.tail=FALSE)
	if(any(pval <= .Machine$double.eps, na.rm=TRUE))
		warning("Some of the p-values are very small (< 2.2e-16). These results might be\n",
			"misleading, as the small p-values might be caused by non-biological artefacts\n",
			"(e.g., sparseness of the data).", call.=FALSE)
	if(is.null(rn1))
		rn1 <- paste("SNP", 1:n.snp1, sep="")
	if(n.snp2 > 1){
		if(is.null(rn2))
			rn2 <- paste("SNP", n.snp1 + (1:n.snp2), sep="")
		names(ll.full) <- names(stat) <- names(pval) <- paste(rep(rn1, e=n.snp2), rep(rn2, n.snp1), sep=" : ")
	}
	else
		names(ll.full) <- names(stat) <- names(pval) <- rn1
	if(!is.null(valMAF)){
		mat.maf <- cbind(rep(round(valMAF[[1]], 4), e=n.snp2), rep(round(valMAF[[2]], 4), n.snp1))
		rownames(mat.maf) <- names(ll.full)
		colnames(mat.maf) <- c("First MAF", "Second MAF")
	}
	else
		mat.maf <- NULL
	out <- list(ll.main=ll.main, ll.full=ll.full, stat=stat, pval=pval,
		maf=valMAF, matMAF=mat.maf)
	class(out) <- "colTDTepi"
	out
}
	
allBetweenCombs <- function(gene){
	ids <- 1:length(gene)
	mat.comb <- NULL
	for(i in 1:(max(gene) - 1)){
		tmp1 <- ids[gene == i]
		tmp2 <- ids[gene > i]
		tmpmat <- cbind(rep(tmp1, e=length(tmp2)), rep.int(tmp2, length(tmp1)))
		mat.comb <- rbind(mat.comb, tmpmat)
	}
	mat.comb
}


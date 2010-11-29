colGxE <- function(mat.snp, env, model=c("additive", "dominant", "recessive"), size=50){
	if (!is.matrix(mat.snp)) 
        	stop("mat.snp has to be a matrix.")
    	if (nrow(mat.snp)%%3 != 0) 
        	stop("mat.snp does not seem to contain trio data, as its number of rows is\n", 
            		"not dividable by 3.")
    	if (is.null(rownames(mat.snp))) 
        	stop("mat.snp does not seem to be a matrix in genotype format,\n", 
            		"as the names of the rows are missing.")
    	if (any(!mat.snp %in% c(0, 1, 2, NA))) 
        	stop("The values in mat.snp must be 0, 1, and 2.")
	if(size<1)
		stop("size should be at least 1.")
	if(length(env)!=nrow(mat.snp)/3)
		stop("The length of env must be equal to the number of trios in mat.snp.")
	if(any(is.na(env)))
		stop("No missing values are allowed in env.")
	if(any(!env %in% 0:1))
		stop("The values in env must be 0 and 1.")
	if(sum(env)<5 | sum(1-env)<5)
		stop("Each of the two groups specified by env must contain at least 5 trios.") 
	type <- match.arg(model)
	fun <- match.fun(switch(type, additive=fastGxEsplit, dominant=fastGxEdomSplit,
		recessive=fastGxErecSplit))
	env2 <- rep(env, e=3)
	tmp1 <- fun(mat.snp, env2==0, size=size)
	tmp2 <- fun(mat.snp, env2==1, size=size)
	beta <- cbind(SNP=tmp1[,1], GxE=tmp2[,1]-tmp1[,1])
	se <- cbind(SNP=sqrt(tmp1[,2]), GxE=sqrt(tmp2[,2]+tmp1[,2]))
	rownames(beta) <- rownames(se) <- colnames(mat.snp)
	stat <- beta/se
	stat <- stat*stat
	lower <- exp(beta - qnorm(0.975) * se)
	upper <- exp(beta + qnorm(0.975) * se)
	pval <- pchisq(stat, 1, lower.tail=FALSE)
	out <- list(coef=beta, se=se, stat=stat, pval=pval, OR=exp(beta), lowerOR=lower, upperOR=upper,
		type=type)
	class(out) <- "colGxE"
	out
}


fastGxEsplit <- function(geno, env2, size=50){
	n.snp <- ncol(geno)
	int <- unique(c(seq.int(1, n.snp, size), n.snp+1))	
	num <- denom <- numeric(n.snp)
	for(i in 1:(length(int)-1)){
		tmp <- fastTDTchunk(geno[env2,int[i]:(int[i+1]-1), drop=FALSE])
		num[int[i]:(int[i+1]-1)] <- tmp$num
		denom[int[i]:(int[i+1]-1)] <- tmp$denom
	}
	logit <- function(x) log(x/(1-x))
	beta <- logit(num/denom)
	se <- denom / ((denom-num)*num)
	cbind(beta, se)
}
	


fastGxEdomSplit <- function(geno, env2, size=50){
	n.snp <- ncol(geno)
	int <- unique(c(seq.int(1, n.snp, size), n.snp+1))
	dmat <- matrix(0, n.snp, 4)
	for(i in 1:(length(int)-1))
		dmat[int[i]:(int[i+1]-1),] <- fastTDTdomChunk(geno[env2,int[i]:(int[i+1]-1), drop=FALSE])
	rownames(dmat) <- colnames(geno)
	h <- (dmat[,1]/3 - dmat[,2] + dmat[,3] - dmat[,4]/3) / (2*(dmat[,1]+dmat[,3]))
	tmp <- (dmat[,2]+dmat[,4])/(3*(dmat[,1]+dmat[,3])) + h*h
	or <- sqrt(tmp) - h
	beta <- log(or)
	tmp <- (dmat[,2]+dmat[,1])*or/(or+1)^2 + (dmat[,3]+dmat[,4])*or/(3*(or+1/3)^2)
	cbind(beta, 1/tmp)
}

fastGxErecSplit <- function(geno, env2, size=50){
	n.snp <- ncol(geno)
	int <- unique(c(seq.int(1, n.snp, size), n.snp+1))
	rmat <- matrix(0, n.snp, 4)
	for(i in 1:(length(int)-1))
		rmat[int[i]:(int[i+1]-1),] <- fastTDTrecChunk(geno[env2,int[i]:(int[i+1]-1), drop=FALSE])
	rownames(rmat) <- colnames(geno)
	h <- (3*rmat[,1] - rmat[,2] + rmat[,3] - 3*rmat[,4]) / (2 * (rmat[,1]+rmat[,3]))
	tmp <- 3 * (rmat[,2] + rmat[,4]) / (rmat[,1] + rmat[,3]) + h*h
	or <- sqrt(tmp) - h		 
	beta <- log(or)
	tmp <- (rmat[,1] + rmat[,2]) * or/(or+1)^2 + 3 *(rmat[,3] + rmat[,4]) * or / (or+3)^2
	cbind(beta,1/tmp)
}


print.colGxE <- function(x, top=5, digits=4, onlyGxE=FALSE, ...){
	if(!onlyGxE){
		pvalG <- format.pval(x$pval[,1], digits=digits)
		outG <- data.frame(Coef=x$coef[,1], OR=x$OR[,1], Lower=x$lowerOR[,1], Upper=x$upperOR[,1],
			SE=x$se[,1], Statistic=x$stat[,1], "p-value"=pvalG, check.names=FALSE)
	}
	pvalGE <- format.pval(x$pval[,2], digits=digits)
	outGE <- data.frame(Coef=x$coef[,2], OR=x$OR[,2], Lower=x$lowerOR[,2], Upper=x$upperOR[,2],
		SE=x$se[,2], Statistic=x$stat[,2], "p-value"=pvalGE, check.names=FALSE)
	cat("          Genotypic TDT for GxE Interactions with Binary E\n\n", "Model Type: ", 
		switch(x$type, "additive"="Additive", "dominant"="Dominant","recessive"="Recessive"), 
		"\n\n",sep="")
	if(nrow(x$coef) > top){
		ord <- order(x$pval[,2])[1:top]
		if(!onlyGxE)
			outG <- outG[ord,]
		outGE <- outGE[ord,]
		cat("Top", top, "GxE Interactions:\n")
	}
	else
		cat("Effects of the GxE Interactions:\n")
	print(format(outGE, digits=digits))
	if(!onlyGxE){
		cat("\n\n", "Effects of the SNPs in the Corresponding GxE Models:\n", sep="")
		print(format(outG, digits=digits))
	}
}

	






	


	
trioLR <- function(x, ...) UseMethod("trioLR")

trioLR.formula <- function(formula, data, recdom=TRUE, ...){
	require(logicFS, quietly=TRUE) || stop("Package logicFS is required.")
	xy <- getXy(formula, data, recdom=recdom)
	trioLR(xy$x, xy$y, ...)
}

trioLR.trioPrepare <- function(x, ...){
	trioLR(x$bin[,-1], x$bin[,1], ...)
}

trioLR.default <- function(x, y, search=c("sa", "greedy", "mcmc"), nleaves=5, penalty=0, weights=NULL,
		control=lrControl(), rand=NA, ...){
	require(LogicReg, quietly=TRUE) || stop("Package LogicReg is required.")
	call <- match.call()
  	if(nrow(x)%%4 != 0)
    		stop("x does not seem to contain trio data with three pseudo-controls for each affected children\n",
      			"since the number of rows is not a multiple of 4.", call.=FALSE)
  	n.trios <- nrow(x)/4
	if(missing(y))
		y <- rep.int(c(3,0,0,0), n.trios)
	if(is.null(weights))
		weights <- rep.int(1, n.trios)
	choice <- match.arg(search)
	select <- checkTLRinput(x, y, choice, n.trios, nleaves=nleaves, penalty=penalty, weights=weights,
		control=control, rand=rand)
	out <- triologreg(x, y, weights, select, nleaves=nleaves, penalty=penalty, control=control, rand=rand)
	out$call <- call
	class(out) <- "trioLR"
	out
}



triologreg <- function(bin, resp, weights, choice, nleaves=5, penalty=0, control=NULL, rand=NA, notBuiltStop=TRUE){
	logreg.storetree <- function(x){
        	i1 <- matrix(x[-(1:3)], ncol = 4, byrow = TRUE)
        	i2 <- length(i1[, 1])
        	i3 <- data.frame(1:i2, i1)
        	names(i3) <- c("number", "conc", "knot", "neg", "pick")
        	i3
    	}
	wgt <- rep(weights, e=4)
	seed <- ifelse(is.na(rand), 0, rand)
	binnames <- colnames(bin)
	mdl <- 9
	n1 <- length(resp)
    	n2 <- ncol(bin)
    	nsep <- 0
	sep <- 0
	nrep <- 25
    	tree.control <- list(treesize=2^floor(log2(control$treesize)), opers=control$opers, 
		minmass=control$minmass)
  	cens <- rep(1, n1)
    	if(choice != 2){
		ntr <- c(1, 1)
        	msz <- rep(nleaves, 2)
	}
    	else{
		ntr <- c(1, 1)
		msz <- nleaves
	}
    	anneal.control <- list(start=control$start, end=control$end, iter=control$iter, 
		earlyout=control$earlyout, update=control$update)
	niter <- ifelse(control$iter==0, 50000, control$iter)
	mc.control <- list(nburn=control$nburn, niter=niter, hyperpars=c(control$hyperpars, 0, 0, 0),
		update=control$update, output=control$output)
	if(choice==7)
		anneal.control$update <- mc.control$update
	if(!is.na(rand))
		xseed <- rand
	else
        	xseed <- floor(runif(1) * 1e+06) + 1
	ipars <- c(mdl, msz[1:2], tree.control$treesize, ntr[1:2], tree.control$opers, 
		anneal.control$update, xseed, nrep, choice, anneal.control$earlyout, 
		tree.control$minmass, mc.control$nburn, mc.control$niter, mc.control$output)
	rpars <- c(anneal.control$start, anneal.control$end, anneal.control$iter, 
		penalty, mc.control$hyperpars)
	orders <- order(rank(resp) + runif(n1)/100)
	nkn <- tree.control$treesize * 2 - 1
    	nxx <- 2
    	if(choice == 1){
		n.a <- ntr[1] * (nkn * 4 + 3)
        	n.b <- nsep + ntr[1] + 1
        	n.c <- 2
    	}
    	if(choice == 2){
        	n.d <- (ntr[2] - ntr[1] + 1) * (msz[2] - msz[1] + 1)
        	n.a <- n.d * ntr[2] * (nkn * 4 + 3)
		n.b <- n.d * (nsep + ntr[2] + 1)
        	n.c <- n.d
    	}
    	if(choice == 6){
        	n.d <- msz[2] + 2
		n.a <- n.d * ntr[2] * (nkn * 4 + 3)
        	n.b <- n.d * (nsep + ntr[2] + 1)
        	n.c <- n.d
    	}
	n100 <- -100
	if(choice == 7){
		n.a <- 256
        	n.b <- n2
        	n.c <- n2 * n2
        	if(abs(mc.control$output) < 2) 
            		n.c <- 1
        	nxx <- n.c * n2
        	if(abs(mc.control$output) < 3) 
            		nxx <- 1
		n100 <- 0
    	}
	xtree <-  rep(n100, n.a)
	ip4 <- 2 * ipars[4] + 1
	fit <- .Fortran("slogreg", as.integer(n1), as.integer(n2), 
		as.integer(nsep), ip = as.integer(ipars), as.single(rpars), 
		as.single(t(sep)), as.integer(cens), as.integer(orders), 
		as.single(resp), as.single(wgt), as.integer(t(bin)), 
		trees = as.integer(xtree), coef = as.single(rep(n100, n.b)), 
		scores = as.single(rep(n100, n.c)), as.integer(ipars[6]), 
		as.integer(ip4), as.integer(rep(0, 2 * ipars[6] * ip4 * n1)), 
		as.integer(rep(0, 7 * ipars[6] * (ip4 + 1) * n2 * 4)), 
		as.single(rep(0, 7 * ipars[6] * (ip4 + 1) * n2 * 4)), 
		as.integer(t(bin)),rd4 = as.integer(rep(0,nxx)),
		PACKAGE = "LogicReg")
	type <- "Trio Logic Regression"
	if(choice %in% (1:2) && is.na(fit$coef[2])){
		if(notBuiltStop)
			cat("For some unknown reason, the fitting of the trio logic regression model failed.\n",
				"Please rerun the analysis with trioLR (maybe after restarting R).\n",
				"If this will not work, please contact the maintainer of the trio package.\n\n", sep="")
		else
			return(NULL)
	}
    	if(choice == 1) 
        	chs <- "single.model"
    	if(choice == 2) 
        	chs <- "multiple.models"
    	if(choice == 7) 
        	chs <- "Bayesian"
    	if(choice == 6) 
        	chs <- "greedy"
	tree.control$operators <- switch(tree.control$opers, "both", "and", "or", "both")
	m1 <- list(nsample = n1, nbinary = n2, nseparate = nsep, type = type, select = chs, 
		seed = rand, choice = choice, control=control)
    	if(choice == 7){
		v1 <- fit$trees
        	v3 <- 1:length(v1)
        	v3 <- max(v3[v1 > 0])
        	v1 <- v1[1:v3]
        	v1 <- data.frame(0:(v3 - 1), v1)
        	names(v1) <- c("size", "number")
        	m1$size <- v1
        	m1$single <- fit$coef
        	m1$single[m1$single < 0] <- 0
        	m1$mc.control <- mc.control
        	if(abs(mc.control$output) > 1){
            		m1$double <- matrix(fit$scores, ncol = length(fit$coef))
            		m1$double[m1$double < 0] <- 0
        	}
        	if(abs(mc.control$output) > 2){
            		m1$triple <- array(fit$rd4, dim = c(length(fit$coef), 
                		length(fit$coef), length(fit$coef)))
		}
	}
    	if(choice == 1){
        	m1$nleaves <- msz[1]
        	m1$ntrees <- ntr[1]
    	}	
    	else{
        	m1$nleaves <- msz
        	m1$ntrees <- ntr
    	}
    	m1$response <- resp
    	m1$binary <- bin
    	m1$separate <- sep
    	m1$censor <- cens
    	m1$weight <- wgt
	if(choice == 1){
		m1$penalty <- penalty
        	m2 <- list()
        	m2$ntrees <- ntr
        	m2$nleaves <- msz
        	m2$score <- fit$scores[1]
        	m2$coef <- fit$coef[1:(nsep + ntr[1] + 1)]
        	class(m2) <- "logregmodel"
        	lx <- 3 + 4 * nkn
        	m2$trees <- list()
        	for(i in 1:ntr[1]){
            		m3 <- list()
            		m3$whichtree <- i
            		m3$coef <- fit$coef[1 + nsep + i]
            		m3$trees <- logreg.storetree(fit$trees[(i - 1) * lx + (1:lx)])
            		class(m3) <- "logregtree"
            		m2$trees[[i]] <- m3
        	}
        	m1$model <- m2
    	}
    	if(choice == 6){
        	m1$nmodels <- fit$ip[1]
        	m1$alltrees <- list()
        	m1$allscores <- matrix(0, nrow = m1$nmodels, ncol = 3)
        	j <- 0
        	k <- 0
        	for(i in 1:m1$nmodels){
            		m2 <- list()
            		m2$score <- fit$scores[i]
            		m2$nleaves <- fit$trees[j + 1]
            		m2$ntrees <- fit$trees[j + 2]
            		m1$allscores[i, ] <- c(fit$scores[i], fit$trees[j + (1:2)])
            		m2$coef <- fit$coef[k + (1:(m2$ntrees + nsep + 1))]
            		class(m2) <- "logregmodel"
            		m2$trees <- list()
            		lx <- 3 + 4 * nkn
            		for(l in 1:m2$ntrees){
                		m3 <- list()
                		m3$whichtree <- l
                		m3$coef <- m2$coef[1 + nsep + l]
                		m3$trees <- logreg.storetree(fit$trees[j + (l - 1) * lx + (1:lx)])
                		class(m3) <- "logregtree"
                		m2$trees[[l]] <- m3
            		}
            		j <- j + lx * m2$ntrees
            		k <- k + (m2$ntrees + nsep + 1)
            		m1$alltrees[[i]] <- m2
        	}
    	}
    	if(choice == 2){
        	m1$nmodels <- fit$ip[1]
        	m1$alltrees <- list()
        	m1$allscores <- matrix(0, nrow = m1$nmodels, ncol = 3)
        	j <- 0
        	k <- 0
        	for(i in 1:m1$nmodels){
            		m2 <- list()
            		m2$score <- fit$scores[i]
            		m2$nleaves <- fit$trees[j + 1]
            		m2$ntrees <- fit$trees[j + 2]
            		m1$allscores[i, ] <- c(fit$scores[i], fit$trees[j + (1:2)])
            		m2$coef <- fit$coef[k + (1:(m2$ntrees + nsep + 1))]
            		class(m2) <- "logregmodel"
            		m2$trees <- list()
            		lx <- 3 + 4 * nkn
            		for(l in 1:m2$ntrees){
                		m3 <- list()
                		m3$whichtree <- l
                		m3$coef <- m2$coef[1 + nsep + l]
                		m3$trees <- logreg.storetree(fit$trees[j + (l - 1) * lx + (1:lx)])
                		class(m3) <- "logregtree"
                		m2$trees[[l]] <- m3
            		}
            		j <- j + lx * m2$ntrees
            		k <- k + (m2$ntrees + nsep + 1)
            		m1$alltrees[[i]] <- m2
        	}
    	}
    	m1$binnames <- binnames
    	m1
}

print.trioLR <- function(x, asDNF=FALSE, posBeta=FALSE, digits=3, ...){
	if(asDNF)
		require(mcbiopi, quietly=TRUE) || stop("Package mcbiopi is required.")
	cat("         Trio Logic Regression\n\n", sep="")
	if(x$choice==7){
		cat("Search Algorithm: MCMC\n\n", sep="")
		if(x$mc.control$output>0)
			cat("Visited models are stored in triolrlisting.tmp.\n")
		else
			cat("Visited models are not stored.\n")
	}
	if(x$choice %in% 1:2){
		cat("Search Algorithm: Simulated Annealing\n\n", sep="")
		if(x$choice==1)
			cat("A single model has been fitted:\n", sep="")
		if(x$choice==2)
			cat(x$nmodels, " models have been fitted:\n", sep="")
	}
	if(x$choice == 6)
		cat("Search Algorithm: Greedy\n\n", x$nmodels, " models have been fitted:\n", sep="")
	cn <- if(is.null(x$binnames)) paste("S", 1:ncol(x$binary), sep="") else x$binnames
	if(x$choice==1){
		printTrioTree(x$model$tree[[1]], cn, asDNF=asDNF, posBeta=posBeta, digits=digits)
		cat("Score: ", signif(x$model$score, digits), "\n")
	}
	if(x$choice %in% c(2,6)){
		score <- signif(x$allscores[,1], digits)
		nleaves <- x$allscores[,2]
		for(i in 1:x$nmodels){
			cat("\n", "Model with ", ifelse(nleaves[i]>0, "a maximum of ", ""), nleaves[i], 
				ifelse(nleaves[i]==1, " leaf", " leaves"), ":\n", sep="")
			printTrioTree(x$alltrees[[i]]$trees[[1]], cn, asDNF=asDNF, posBeta=posBeta, digits=digits)
			cat("Score: ", score[i], "\n", sep="")
		}
	}
}

printTrioTree <- function(ltree, cn, asDNF=FALSE, posBeta=FALSE, digits=3){
	mat <- as.matrix(ltree$trees)
	coef <- signif(ltree$coef, digits)
	n.row <- nrow(mat)
	vec <- character(nrow(mat))
	ids3 <- which(mat[,2] == 3)
	if(length(ids3)==0)
		cat("Model contains no variables.\n", sep="")
	if(length(ids3)==1)
		cat(coef, " * ", ifelse(mat[1,4]==1, "!", ""), cn[mat[1,3]], "\n", sep="")
	if(length(ids3)>1){
		if(asDNF){
			if(posBeta)
				coef <- abs(coef)
			else
				ltree$coef <- abs(ltree$coef)
			tmppi <- getPImps(ltree, 1)
			tmpdnf <- paste("(", getSNPdummyNames(tmppi, cn), ")", sep="", collapse=" | ")
			vec[1] <- if(length(tmppi)==1) tmpdnf else paste("(", tmpdnf, ")", sep="")
		}
		else{			
			for(i in ids3)
				vec[i] <- paste(ifelse(mat[i,4]==1, "!", ""), cn[mat[i,3]], sep="")
			nhalf <- floor(n.row/2)
			for(i in nhalf:1){
				if(mat[i,2] == 1)
					vec[i] <- paste("(", vec[2*i], " & ", vec[2*i+1], ")", sep="")
				if(mat[i,2] == 2)
					vec[i] <- paste("(", vec[2*i], " | ", vec[2*i+1], ")", sep="")
			}
		}
		cat(coef, " * ", vec[1], "\n", sep="")
	}
}


plot.trioLR <- function(x, whichTree=NA, freqType=1, useNames=FALSE, addStats=TRUE, digits=3, main=NULL, cexOper=1.5, 
		cexLeaf=1.5, sizeLeaf=7, cexPar=1.3, ...){
	if(useNames){
		if(is.null(x$binnames))
			stop("If useNames = TRUE, the column names of the genotype matrix must have been specified.")
		cn <- x$binnames
	}
	else
		cn <- 1:ncol(x$binary)
	if(x$choice==1)
		plotTrioTree(x$model$trees[[1]], cn, x$model$score, addStats=addStats, digits=digits, main=main, cexOper=cexOper,
			cexLeaf=cexLeaf, sizeLeaf=sizeLeaf, cexPar=cexPar)
	if(x$choice %in% c(2,6)){
		if(is.na(whichTree))
			stop("If multiple models have been fitted, whichTree must be specified.")
		if(!whichTree %in% 1:x$nmodels)
			stop("whichTree must be an integer between 1 and ", x$nmodels, ".")
		plotTrioTree(x$alltrees[[whichTree]]$trees[[1]], cn, x$allscores[whichTree, 1], addStats=addStats,
			digits=digits, main=main, cexOper=cexOper, cexLeaf=cexLeaf, sizeLeaf=sizeLeaf, cexPar=cexPar)
	}
	if(x$choice == 7){
		if(!freqType %in% 1:3)
			stop("For plotting the results of a MC Trio Logic Regression, freqType must be set\n",
				"either to 1 (for plotting the percentage of visited models containing the individual variables),\n",
				"or to 2 (for plotting the percentage of models containing the different pairs of variables),\n",
				"or to 3 (for the observed-to-expected ratio for the pairs being jointly in the models).")
		if(freqType==1){
			if(length(x$single) == 0)
				stop("No information available in x on the number of models containing the binary variables.")
			if(is.null(main))
				main <- "Marginal Frequency of Being in the Models"
			las <- ifelse(useNames, 3, 0)
			xlab <- ifelse(useNames, "", "Variable")
			barplot(as.vector(x$single)/sum(x$size[,2]) * 100, names=cn, main=main, xlab=xlab,
				ylab="Percentage", las=las)
		}
		else{
			if(length(x$double) == 0)
				stop("No information available in x on the number of models containing the pairs of variables.")
			if(is.null(main))
				main <- paste(ifelse(freqType==2, "Frequency", "Observed/Expected Ratio"), "of Being Jointly in the Models")
			tmp2 <- x$double
			tmp3 <- 1:length(tmp2[,1])
			tmp4 <- outer(tmp3, tmp3, "-")
			tmp2[tmp4 <= 0] <- NA
			if(freqType==2)
				image(tmp3, tmp3, tmp2, main=main, xlab="", ylab="",xaxt="n", yaxt="n")
			else{
				tmp <- x$single
				tmp5 <- outer(tmp/sum(tmp), tmp)
				image(tmp3, tmp3, tmp2/tmp5, main=main, xlab="", ylab="", xaxt="n", yaxt="n")
			}
			axis(1, at=1:length(tmp3), labels=cn, las=ifelse(useNames, 2, 0))
			axis(2, at=1:length(tmp3), labels=cn, las=ifelse(useNames, 1, 0))
		}
	}
}

plotTrioTree <- function(ltree, cn, score, addStats=TRUE, digits=3, main=NULL, cexOper=1.5, cexLeaf=1.5, 
		sizeLeaf=7, cexPar=1.3){
	dat <- ltree$trees
	coef <- signif(ltree$coef, digits)
	score <- signif(score, digits)
	if(is.null(main))
		main <- "Logic Tree"
	if(dat[1,5] == 0){
		cat("The trio logic regression model contains no variable.\n")
		plot(0:1, 0:1, type="n", xaxt="n", yaxt="n", xlab="", ylab="", main=main)
		text(0.5, 0.5, "Empty Tree", cex=2)
		if(addStats)
			title(sub=paste("Score:", score, "\n"), cex.sub=cexPar)
	}
	else{
		level <- ceiling(log2(1 + max(which(dat[,5]==1))))
		max.p <- 2^level - 1
		dat <- dat[1:max.p,]
		conc <- dat[,2]
		cn <- c("", cn)
		knot <- cn[dat[,3]+1]
		neg <- dat[,4]
		num <- dat[,1]
		pick <- dat[,5]
		x <- y <- rep(0, max.p)
		left.b <- 2^(0:(level-1))
		right.b <- 2^(1:level)-1
		l.max <- max(num[left.b][pick[left.b]==1])
		r.max <- max(num[right.b][pick[right.b]==1])
		l.shift <- 1/(l.max+1)
		r.shift <- 1/(r.max+1)
		w.shift <- min(l.shift, r.shift)
		plot(c(0.85 * l.shift, 0.95*(1-r.shift)), c(-1, 1.05*(-level)), type="n", xlab="", ylab="", xaxt="n", yaxt="n", col=0, main=main)
		if(addStats)
			title(sub = paste("Parameter:", coef, "      Score:", score, "\n"), cex.sub=cexPar)
		for(k in level:1){
			y.val <- (-k)
			min.val <- 2^(k-1)
			max.val <- 2^k - 1
			diff.x <- 1/(min.val + 1)
			x.val <- seq(diff.x, 1-diff.x, diff.x)
			if(level==1)
				x.val <- x.val * (0.95 * (1-r.shift) - 0.85 * l.shift) + 0.85 * l.shift
			for(j in 1:length(x.val)){
				ids <- 2^(k-1) + j - 1
				if(pick[ids] == 1){
					x[ids] <- x.val[j]
					y[ids] <- y.val
					if(conc[ids] < 3)
						text(x.val[j], y.val, ifelse(conc[ids]==1, "AND", "OR"), cex=cexOper, font=2)
					else{
						points(x.val[j], y.val, pch=ifelse(neg[ids]==0, 0, 15), cex=sizeLeaf)
						text(x.val[j], y.val, knot[ids], cex=cexLeaf, col=ifelse(neg[ids]==0, 1, "white"))
					}
				}
			}
		}
		for(k in 1:max.p){		
			if((ceiling(log2(k)) < level) && (2*k<max.p) && (pick[k] + pick[2*k] == 2)){
				lines(x[c(k, 2*k)], c(y[k] - 0.2, y[2*k] + 0.2))
				lines(x[c(k, 2*k+1)], c(y[k]-0.2, y[2*k+1] + 0.2))
			}
		} 
		invisible()
	}
}				


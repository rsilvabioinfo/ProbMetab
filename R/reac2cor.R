#' reac2cor 
#'
#' Use the intensity of putative molecules in repeated samples 
#' to calculate correlations and partial correlation in a user
#' defined threshold of false discovery rate for significance testing.
#' After the correlation test the function also overlay significant correlations
#' with all putative reactions between two masses.
#'
#' @param mw two column of adjacency matrix indexes connecting compounds
#' by reactions.
#' @param classTable classification table, with intensities for repeated samples.  
#' @param opt correlation option, "cor" for correlation,
#' and "pcor" for partial correlation. 
#' @param corths correlation intensity threshold. 
#' @param corprob probability that the correlation is considered significant. 
#' @param pcorprob probability that the partial correlation is considered significant. 
#'  
#' @return A list of estimated correlations and reactions. 
#'
#' @export


reac2cor <- function(mw, classTable, opt="cor", corths=0.75, corprob=0.8, pcorprob=0.8) {
#	require(GeneNet)

	# for matrices with the same number of columns
	comp.matrix <- function (m1, m2) {
		idx <- duplicated(rbind(m1, m2))
		return(rbind(m1, m2)[idx,])
	}
	
	# very fast way to do the comparation, need to understand the indexes
#	comp.matrix2 <- function(m1, m2) {	
#		require(data.table)
#		M1 = setkey(data.table(m1))
#		M2 = setkey(data.table(m2))
#		idx <- na.omit(
#    			M2[M1,which=TRUE]
#		)
#		return(m1[idx,])
#	}
#

	int.matrix <- classTable[classTable[,6]!="", 8:ncol(classTable)]
	#int.matrix <- int.matrix[-1,] 
	int.matrix <- t(int.matrix) 
	int.matrix <- apply(int.matrix, 2, as.numeric) 
	colnames(int.matrix) <- paste("C", sprintf("%004d", 1:ncol(int.matrix)), sep="")
	rownames(int.matrix) <- as.matrix(classTable[1,][8:ncol(classTable)]) 
	fac <- as.vector(classTable[,1]) 
	fac[fac!=""] <- 1:length(fac[fac!=""]) 
	for(i in 1:length(fac)) if(fac[i]=="") fac[i] <- fac[i-1]
	fac <- cbind(fac, 1:length(fac)) 
	fac <- apply(fac, 2, as.numeric) 
	mwb <- mw 
	mw <- as.matrix(mw) 
	idx <- mw[,1] < mw[,2]
	mw <- mw[idx,] 
	# as a mass (putative metabolite) has many candidates, the set of compound candidates
	# has to map to a set of unique masses
	apply(fac, 1, function(x) { if(length(which(mw==x[2]))) mw[which(mw==x[2], arr.ind=TRUE)] <<- x[1] })
	mw <- apply(mw, 2, as.numeric) 
	# following GeneNet output the undirected nodes are reported from the node of small index to the node of bigger index
	mw2 <- unique(mw) 

	if(opt=="pcor") {
		inferred.pcor <- ggm.estimate.pcor(int.matrix) 
		test.results <- ggm.test.edges(inferred.pcor, plot=FALSE)
		# user set threshold of false discovery rate for nodes
		signif <- (test.results$prob > pcorprob) 
		test.results2 <- test.results[signif,] 

		pcor.vs.reac <- comp.matrix(mw2, as.matrix(test.results2[,2:3]))
		return(list(int.matrix=int.matrix, signif.pcor=test.results2, node.reacs=mw2, 
			pcor.vs.reac=pcor.vs.reac
		))
	}


	# calculation of cor, with significance testing and
	# false discovery rate for multiple testing

	if(opt=="cor") {
		cor.prob <- function(X, dfr = nrow(X) - 2) { 
		  R <- cor(X) 
		  above <- row(R) < col(R) 
		  r2 <- R[above]^2
		  Fstat <- r2 * dfr / (1 - r2)
		  R[above] <- 1 - pf(Fstat, 1, dfr) 
		  R
		} 
		pcor1 <- cor.prob(int.matrix) 
		above <- row(pcor1) < col(pcor1)
		test.results5 <- cbind(t(pcor1)[above], which(above, arr.ind=TRUE), pcor1[above])
		fdr <- fdrtool(test.results5[,4], statistic="pvalue",plot=FALSE)
		test.results5 <- cbind(test.results5, fdr$qval, 1-fdr$qval) 
		cnames <- c("cor", "node1", "node2", "pval", "qval", "prob")
		colnames(test.results5) <- cnames 
		signif <- which(abs(test.results5[,1]) > corths & test.results5[,6]>  corprob)
		test.results5 <- test.results5[signif,]
		 
		# reference mapping of rownames from classTable to crossreference with
		# correlation output
		cor.vs.reac <- comp.matrix(mw2, test.results5[,2:3])  

	#	mw6 <- comp.matrix(cor.vs.reac, pcor.vs.reac) 
	#	mw7 <- comp.matrix(test.results2[,2:3], test.results5[,2:3]) 

		return(list(int.matrix=int.matrix, signif.cor=test.results5, node.reacs=mw2, 
			cor.vs.reac=cor.vs.reac  
		))
	}

}

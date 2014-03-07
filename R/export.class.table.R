#' export.class.table
#'
#' Builds a matrix with the probability for all mass to candidate compounds
#' assignments, by averaging the number of assignments obtained by the gibbs sampler algorithm
#' or ordering the compound candidates with the likelihood matrix. 
#'
#' @param gibbsL a list of attributions and probabilities from gibbs.samp function. 
#' @param reactionM data.frame with compound annotation information. 
#' @param molIon non redundant ion annotation. 
#' @param probM optionally to gibbsL, a matrix of likelihoods. 
#' @param html logical, indicating whether a html file should be generated. This parameter uses the raw data to plot EICs and may be time consuming.
#' @param filename html file name, the default is "test". 
#' @param burnIn how many samples of the gibbs sampler should be discarded.
#' @param linkPattern which pattern should be linked to compound id, for now we have 
#' implemented "kegg", "pubchem" and "chebi" patterns.
#' @param m.test statistical test to compare mean differences. This option
#' is only available to single acquisition mode analysis, with options
#' "t.test" and "anova". 
#' @param class1 if the m.test is "t.test" first class to compare in the test,
#' according with xcmsSet phenoData.
#' @param class2 if the m.test is "t.test" second class to compare in the test,
#' according with xcmsSet phenoData.
#' @param norm logical, if TRUE performs median normalization from (Anal. Chem. 2011, 83, 5864-5872).
#' @param DB data.frame table used to search compounds, with the field name to be incorporated in the html table. 
#' @param prob how to calculate the probability to attribute a mass to a compound. 
#' Default is "count", which divide the number of times each identity was 
#' was attributed by the number of samples. Optionally the user could
#' choose to use the mean of the probabilities of the identity, "mean". 
#' @return A list with a matrix "classTable" with attributions and probabilities and
#' indexes of selected masses from xcms peak table. 
#'
#' @export

export.class.table <- function(gibbsL=NULL, reactionM, molIon=NULL, probM=NULL, html=FALSE, filename="test", burnIn=3000, linkPattern="kegg", m.test="none", class1=NULL, class2=NULL, norm=FALSE, DB, prob="count") {       

	plotEIC <- function (xcmsObject, figidx, pngidx, mode=NULL) {
		dir.create(paste(filename,"_fig",sep=""))
		gt<-groups(xcmsObject)
		gt <- gt[figidx,]
		rgt <- gt[,c("rtmin","rtmax")] 
			rgt[,1] <- rgt[,1]-100 
			rgt[,2] <- rgt[,2]+100 
			eics <- getEIC(xcmsObject, mzrange=gt, rtrange =rgt, groupidx = 1:nrow(gt)) 
			png(file.path(paste(filename, "_fig/%003d.png", sep="")), height=768, width=1024) 
			#png(file.path(paste(filename, "_fig/", pngidx, sep="")), h=768, w=1024) 
			plot(eics, xcmsObject)
			dev.off()
			if(!is.null(mode)) {
				pngs <- dir(paste(filename, "_fig/", sep=""))
				if(length(grep("pos|neg" , pngs))) pngs <- pngs[-grep("pos|neg" , pngs)] 
				opng <- as.numeric(sub(".png","", pngs)) 
				pngs <- pngs[order(opng)]
				name1 <- paste(filename, "_fig/",pngs, sep="") 
				name2 <- paste(filename, "_fig/",pngidx, mode, ".png", sep="")
				for(i in 1:length(name1)) file.rename(name1[i], name2[i])
			}

	}
    allion <- molIon$molIon[molIon$molIon[,"isotope"]==0,] 
    ReactMatrix <- reactionM[reactionM[,5]!="unknown",]
    x <- apply(unique(ReactMatrix[,c(2, 3)]), 2, as.numeric) # Have to look for all pairs
    y <- as.numeric(ReactMatrix[,4])
    prob_mean_ma <- matrix(0, nrow = length(y), ncol = nrow(x))
#    z_average <- matrix(0, nrow = length(y), ncol = length(x))

    if (!is.null(gibbsL)){	
        prob_table <- gibbsL$prob_table[,-c(1:burnIn)]  
        class_table <- gibbsL$class_table[,-c(1:burnIn)]
	#indList <- tapply(1:nrow(ReactMatrix), as.numeric(ReactMatrix[,1]), function(x) x) 
	coords <- tapply(1:nrow(ReactMatrix), ReactMatrix[,"molIonID"], function(x) x)
	coords2 <- unlist(lapply(coords, function(x) rep(x[1], length(x)))) 
	indList <- coords[order(unique(coords2))] 
	fillMatrix <- function(j,i) { 
			idp <- which(class_table[i,] == j)
			if(prob=="count") prob_mean_ma[j,i] <<- length(idp)/ncol(class_table)
			if(prob=="mean") prob_mean_ma[j,i] <<- mean(prob_table[i,idp]) 
		}


        for ( i in 1:nrow(x) ) {
		
			sapply(indList[[i]], fillMatrix, i)
	}
	if(sum(prob_mean_ma=="NaN"))	prob_mean_ma[prob_mean_ma=="NaN"] <- 0
#        for ( i in 1:nrow(x) ) {
#            for ( j in 1:length(y) ) {
#                idp <- which(class_table[i,] == j)
#                prob_mean_ma[j,i] <- mean(prob_table[i,idp]) 
#		# this is an alternative way to calculate the probabilities, should try latter, and compare results
#                #prob_mean_ma[j,i] <- length(idp)/ncol(class_table)
#                if ( prob_mean_ma[j,i] == "NaN" ) prob_mean_ma[j,i] <- 0
#            }
#           # do I still need this matrix? 
#            k=which(prob_mean_ma[,i]==max(prob_mean_ma[,i]))
#            z_average[k[1],i]=1
#        }
    }
    else {
        prob_mean_ma <- probM 
    }
	# think about natural probabilities
	# prob_mean_ma[prob_mean_ma[,1]!=0,1]/sum(prob_mean_ma[prob_mean_ma[,1]!=0,1]) 
	prob_mean_ma <- apply(prob_mean_ma, 2, function(x){ x[x!=0] <- x[x!=0]/sum(x[x!=0]); return(x)} )

	# create a dir to figures
	lpattern <- function(type){
		switch(type,
			kegg = "http://www.genome.jp/dbget-bin/www_bget?",
			chebi = "http://www.ebi.ac.uk/chebi/searchId.do;EFB7DFF9E88306BBCD6AB78B32664A85?chebiId=",
			pubchem = "http://www.ncbi.nlm.nih.gov/pccompound/?term="
		)
	}
        linkURL <- lpattern(linkPattern)
	fig <- paste("file://", getwd(), paste("/",filename,"_fig/",sep=""), sep="")
	if(!is.null(molIon$cameraobj)) {
		figidx <- c("")
		coords <- gsub("(^\\d)","X\\1",rownames(molIon$cameraobj@xcmsSet@phenoData)) 
		# experimental! Which set of characters????
		coords <- gsub("-|\\,|~","\\.",coords)
		coords <- gsub("\\s+","\\.",coords)
		peaklist <- getPeaklist(molIon$cameraobj)
    		rpeaklist <- peaklist[,c("mz","rt","isotopes","adduct","pcgroup")]
	}
	else {
		figidx <- c("","")
		coordsP <- gsub("(^\\d)","X\\1",rownames(molIon$pos@xcmsSet@phenoData)) 
		# experimental! Which set of characters????
		coordsP <- gsub("-|\\,|~","\\.",coordsP)
		coordsP <- gsub("\\s+","\\.",coordsP)
		coordsN <- gsub("(^\\d)","X\\1",rownames(molIon$neg@xcmsSet@phenoData)) 
		# experimental! Which set of characters????
		coordsN <- gsub("-|\\,|~","\\.",coordsN)
		coordsN <- gsub("\\s+","\\.",coordsN)
		coords <- coordsP
		if(length(coordsP)!=length(coordsN)) cat("\n Warning: The number of samples are different\n")
	
		peaklistP <- getPeaklist(molIon$pos)
    		rpeaklistP <- peaklistP[,c("mz","rt","isotopes","adduct","pcgroup")]
		peaklistN <- getPeaklist(molIon$neg)
    		rpeaklistN <- peaklistN[,c("mz","rt","isotopes","adduct","pcgroup")]
	}
	
#	if(sum(is.na(peaklist))) {
#		cat("\nWarning: NAs Found in peaklist\n\nSubstituting for \"ones\"\n")
#		na.ids <- which(is.na(peaklist),arr.ind=TRUE)
#		for(l in 1:nrow(na.ids)){
#			peaklist[na.ids[l,][1], na.ids[l,][2]] <- 1
#		}
#	}
#					

    ans <- matrix("", nrow=1, ncol=7+length(coords))
    unq <- unique(ReactMatrix[,2:3])
        for (i in 1:nrow(unq)) {
		coord <- which(ReactMatrix[,2]==unq[i,1] & ReactMatrix[,3]==unq[i,2])
		coord2 <- which(allion[,2]==unq[i,1] & allion[,1]==unq[i,2])
		#	idx2 <- unique(which(allion[,1] %in% reactionM[reactionM[,5]=="unknown",2]))
		# work with the higher intensities for a given ion annotation, not necessarily the right one
		
		if(!is.null(molIon$cameraobj)) {
		idx <- as.vector(unlist(sapply(allion[coord2,"trace"], 
						function(x) {
							x <- as.matrix(x)
							raw <- strsplit(x,";")[[1]]
							mraw <- apply(peaklist[raw, coords], 1, mean)
							raw[which.max(mraw)]
						}

						)
					)
				) 
		
		idx <- unique(idx) 
		figidx <- append(figidx,idx)
		}
		else {
	idx <- c()

 for(l in 1:nrow( allion[coord2,c("trace","comb")])) { 
						x <- as.matrix(allion[coord2,c("trace","comb")][l,])
							raw <- strsplit(x[1],";")[[1]]
							if(x[2]!="neg"){
							mraw <- apply(peaklistP[raw, coordsP], 1, mean, na.rm=TRUE)
							}
							else {
							
							mraw <- apply(peaklistN[raw, coordsN], 1, mean, na.rm=TRUE)
							}
							idx <-  c(idx, raw[which.max(mraw)])
						}

	
	idx <- unique(idx) 
		figidx <- rbind(figidx,c(idx,allion[coord2,"comb"][1]))
		}
		#figidx <- append(figidx,strsplit(allion[coord2,5], ";")[[1]][1])
    		ans1 <- matrix("", nrow=length(coord), ncol=7+length(coords))
		ans1[,2]<-as.matrix(ReactMatrix[coord,5])
		prob <- as.matrix(prob_mean_ma[coord, i]) # need to change and compare a pair of mass/rt
		# number figs
		if ( i >= 100 ) { ans1[1,6]=i }
		else { if ( i >= 10 ) { ans1[1,6]=paste(0,i, sep="") } else { ans1[1,6]=paste("00",i, sep="") }  } 

		if (sum(prob)>0) {
		    #prob <- prob/sum(prob) 
		    o <- order(prob, decreasing=TRUE)
		    ans1[,-6] <- ans1[o,-6]
		    ans1 <- matrix(ans1, nrow=length(o))
		    ans1[1,1] <- ReactMatrix[coord[1],3]
		    #ans1[,3] <- round(prob/min(prob[prob!=0]), 3)[o] 
		    ans1[,3] <- round(prob, 3)[o] 
		    if (length(prob[prob!=0])>1) { 
			entropy <-  -sum(prob[prob!=0]*log(prob[prob!=0], length(prob[prob!=0]))) 
		    } 
		    else { entropy <- 0 
		    }
		    ans1[1,4] <- round(entropy, 3)
		}
		else {
		    ans1[1,1] <- ReactMatrix[coord[1],3]
		    ans1[1,3] <- "undef"
		}

		if(!is.null(molIon$cameraobj)) {
			ans1[1,7] <- apply(rpeaklist[idx,], 1, function(x) paste(x[c(1,2,3,4)], collapse="#"))
			ans1[1,8:ncol(ans1)] <- as.matrix(peaklist[idx, coords])
		}
		else {
			if(allion[coord2,"comb"]=="pos"|allion[coord2,"comb"]=="both") {
				ans1[1,7] <- apply(rpeaklistP[idx,], 1, function(x) paste(x[c(1,2,3,4)], collapse="#"))
				ans1[1,8:ncol(ans1)] <- as.matrix(peaklistP[idx, coordsP])
			}
			else {
				ans1[1,7] <- apply(rpeaklistN[idx,], 1, function(x) paste(x[c(1,2,3,4)], collapse="#"))
				ans1[1,8:ncol(ans1)] <- as.matrix(peaklistN[idx, coordsN])
			}
		}
		ans <- rbind(ans, as.matrix(ans1))
      }
	ans <- ans[-1,]
	# this option should change according with the bank
	if(html) {
		nid <- unlist(sapply(ans[,2], function(x) which(DB$id==x)))
		#ans[,2] <- as.character(DB$name[nid]) 
	}
	unk <- reactionM[reactionM[,5]=="unknown",]
	ans1 <- matrix("", nrow=nrow(unk), ncol=7+length(coords)) 
	ans1[,1] <- unk[,3]
	ans1[,2] <- unk[,5]
	for(j in 1:nrow(ans1)) {
		i <- j + max(as.numeric(ans[,6]),na.rm=TRUE) 
		if ( i >= 100 ) { ans1[j,6]=i }
		else { if ( i >= 10 ) { ans1[j,6]=paste(0,i, sep="") } else { ans1[j,6]=paste("00",i, sep="") }  } 
	}
	# this step try to recover ids of ion annotation for masses without database annotation
	idx2 <- c(); #for(m in 1:nrow(allion))  if(sum(allion[m,2]==as.numeric(unk[,2])) & sum(allion[m,1]==as.numeric(unk[,3]))) idx2 <- append(idx2, m)
	# temp changes made 03/03/2014 have to check carefuly
	lidx <- lapply(1:nrow(allion), function(m) which(allion[m,2]==unk[,2] & allion[m,1]==unk[,3])) 
	idx2 <- which(lapply(lidx, length)>0) 

	if(!is.null(molIon$cameraobj)) {
	idx <- as.vector(unlist(sapply(allion[idx2,"trace"], 
					function(x) {
						x <- as.matrix(x)
						raw <- strsplit(x,";")[[1]]
						mraw <- apply(peaklist[raw, coords], 1, mean)
						raw[which.max(mraw)]
					}

					)
				)
			) 
	}

	else {
	# don't know what happened here with apply
	idx <- c()

 for(i in 1:nrow( allion[idx2,c("trace","comb")])) { 
						x <- as.matrix(allion[idx2,c("trace","comb")][i,])
							raw <- strsplit(x[1],";")[[1]]
							if(x[2]!="neg"){
							mraw <- apply(peaklistP[raw, coordsP], 1, mean, na.rm=TRUE)
							}
							else {
							
							mraw <- apply(peaklistN[raw, coordsN], 1, mean, na.rm=TRUE)
							}
							idx <-  c(idx, raw[which.max(mraw)])
						}

	

	tmpidx <- cbind(idx,allion[idx2,"comb"])
	}	
	if(!is.null(molIon$cameraobj)) {
		ans1[,7] <- apply(rpeaklist[idx,], 1, function(x) paste(x[c(1,2,3,4)], collapse="#"))
		ans1[,8:ncol(ans1)] <- as.matrix(peaklist[idx, coords])
	}
	else {
		idxP <- tmpidx[tmpidx[,2]!="neg",1]  
		ans1[1:length(idxP),7] <- apply(rpeaklistP[idxP,], 1, function(x) paste(x[c(1,2,3,4)], collapse="#"))
		ans1[1:length(idxP),8:ncol(ans1)] <- as.matrix(peaklistP[idxP, coordsP])
		idxN <- tmpidx[tmpidx[,2]=="neg",1]  
		ans1[(length(idxP)+1):nrow(ans1),7] <- apply(rpeaklistN[idxN,], 1, function(x) paste(x[c(1,2,3,4)], collapse="#"))
		ans1[(length(idxP)+1):nrow(ans1),8:ncol(ans1)] <- as.matrix(peaklistN[idxN, coordsN])
	}
	ans <- rbind(ans, as.matrix(ans1))
	
	if(!is.null(molIon$cameraobj)) {
		figidx <- c(figidx,idx)
		figidx <- as.numeric(figidx[-1])
	}
	else {
		figidx <- rbind(figidx,tmpidx)
		allidx <- figidx[-1,] 
		allidx <- cbind(allidx, ans[ans[,6]!="",6])
		figidx <- as.numeric(figidx[-1,1])
	}

	
	if(m.test=="none") {
		testname <- "none"
		#testname <- "Formula"
		#ans[ans[,2]!="unknown",][,5] <- as.character(DB$formula[nid]) 
	}
	if(m.test=="t.test") {
	normalize.medFC <- function(mat) {
	   # Perform median fold change normalisation
	   #           X - data set [Variables & Samples]
	   medSam <- apply(mat, 1, median)
	   medSam[which(medSam==0)] <- 0.0001
	   mat <- apply(mat, 2, function(mat, medSam){
	      medFDiSmpl <- mat/medSam
	      vec<-mat/median(medFDiSmpl)
	      return(vec)
	   }, medSam)
	   return (mat)
	}
	# this piece of code was copied from xcms
	pval <- function(X, classlabel, teststat) {

	    n1 <- rowSums(!is.na(X[,classlabel == 0]))
	    n2 <- rowSums(!is.na(X[,classlabel == 1]))
	    A <- apply(X[,classlabel == 0], 1, sd, na.rm=TRUE)^2/n1 ## sd(t(X[,classlabel == 0]), na.rm = TRUE)^2/n1
	    B <- apply(X[,classlabel == 1], 1, sd, na.rm=TRUE)^2/n2 ## sd(t(X[,classlabel == 1]), na.rm = TRUE)^2/n2
	    df <- (A+B)^2/(A^2/(n1-1)+B^2/(n2-1))

	    pvalue <- 2 * (1 - pt(abs(teststat), df))
	    invisible(pvalue)
	}

	c1 <- grep(class1, molIon$cameraobj@xcmsSet@phenoData[,1]) 
	c2 <- grep(class2, molIon$cameraobj@xcmsSet@phenoData[,1]) 
	testclab <- c(rep(0, length(c1)), rep(1, length(c2)))
	testval <- groupval(molIon$cameraobj@xcmsSet, "medret", "into") 
	if(norm) testval <- normalize.medFC(testval) 	
	tstat <- mt.teststat(testval, testclab)
        pvalue <- pval(testval, testclab, tstat)

#
#		rport <- diffreport(molIon$cameraobj@xcmsSet, class1=class1, class2= class2, sortpval=FALSE) 
#		ans[ans[,6]!="",5] <- rport[figidx, "pvalue"]
		ans[ans[,6]!="",5] <- pvalue[figidx]
		testname <- "t.test p-value"
	}
	if(m.test=="anova"){
		class <- molIon$cameraobj@xcmsSet@phenoData 
		getPvalue <- function(dataidx) {
			aov.data <- data.frame(resp=as.numeric(peaklist[dataidx,coords]), class=class) 
			anova(aov(resp~class, aov.data))$Pr[1] 
		}
		testname <- "anova p-value"
		ans[ans[,6]!="",5] <- sapply(figidx, getPvalue)
	}

	header <- matrix(c("Proposed Mass","Most probable Compound","Probability","Entropy", testname,"EIC-plot", "Ion annotation",coords), nrow=1 , ncol=7+length(coords) ) 
	ans <- rbind(header, ans)

	# additional field
	# ans <- cbind(ans[,1:2], ans[,2], ans[,3:ncol(ans)]) 
	#ans[ans[,3]!="unknown",][-1,3] <- as.character(DB$sbml.id[nid]) 

	if(html) {
		#require(hwriter)
		ansb <- ans
		ans[ans[,2]!="unknown",][-1,2] <- as.character(DB$name[nid]) 
		if(linkPattern=="pubchem") ansb <- ans
		
		hyper=matrix(paste(linkURL, ansb[-1,2], sep=""),ncol=1 )
		if(!is.null(molIon$cameraobj)) {	
			hyper1=matrix(paste(fig, ans[-1,6],".png", sep=""),ncol=1 ) 	
		}
		else {
			hyper1 <- ans[-1,6]
			hyper1[ans[-1,6]!=""][allidx[,2]!="neg"] <- paste(hyper1[ans[-1,6]!=""][allidx[,2]!="neg"], "pos", sep="")
			hyper1[ans[-1,6]!=""][allidx[,2]=="neg"] <- paste(hyper1[ans[-1,6]!=""][allidx[,2]=="neg"], "neg", sep="")
			hyper1=matrix(paste(fig, hyper1,".png", sep=""),ncol=1 ) 	
		}	
		p=openPage(paste(filename,".html",sep=""))
		ans2 <- ans[,1:7]
		link <- cbind(matrix(NA,nrow(ans2),1),rbind(NA,hyper),matrix(NA,nrow(ans2),3),rbind(NA,hyper1),matrix(NA,nrow(ans2),1))
		hwrite(ans2, p,row.bgcolor='#ffdc98', link=link  )
		closePage(p)
		if(!is.null(molIon$cameraobj)) {	
			plotEIC(molIon$cameraobj@xcmsSet, figidx, ans[ans[,6]!="",6][-1]) 
		}	
		else {
			 dataidxP <- as.numeric(allidx[allidx[,2]!="neg",1]) 
			 pngidxP <- allidx[allidx[,2]!="neg",3]
			plotEIC(molIon$pos@xcmsSet, dataidxP, pngidxP, "pos") 
			 dataidxN <- as.numeric(allidx[allidx[,2]=="neg",1]) 
			 pngidxN <- allidx[allidx[,2]=="neg",3]
			plotEIC(molIon$neg@xcmsSet, dataidxN, pngidxN, "neg") 
		}
		
	}
	else {
		ansb <- ans
	}
	colnames(ansb) <- ansb[1,]
	ansb <- ansb[-1,]
	return(list(classTable=ansb,  figidx=figidx))
} 

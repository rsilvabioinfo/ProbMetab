#' get.annot 
#'
#' This function extracts annotation from CAMERA object,
#' generating a matrix of non-redundant putative molecular ions.
#' @param xsAnnotate CAMERA's annotation object. 
#' @param polarity acquisition mode of mass spectrometer. 
#' @param allowMiss logical, optionally retrieves peaks with no evidence
#' of adduct or isotope and annotate them as single charged molecules [M+/-H]. 
#' @param toexclude samples to be excluded of peak counting to non-annotated peak selection. 
#' @param xset xcmsSet xcms object after missing data replacement,
#' to retrieve SNR to isotopic peaks.
#' @param minsamp minimum number of samples in which an ion should be
#' present to be selected. 
#' @param minint minimum mean intensity that a ion should
#' present to be selected. 
#' @return A list with a matrix of possible molecular ions with a trace of their annotation, and the used xsAnnotate object. 
#'
#' @export

get.annot <- 
function(xsAnnotate, polarity="positive",allowMiss=FALSE, xset=NULL, toexclude=NULL, minsamp=0.6, minint=5000) {

	if(polarity!="positive" & polarity!="negative") stop("Unknown polarity") 
          isoID <- as.matrix(xsAnnotate@isoID)
          annoGrp <- as.matrix(xsAnnotate@annoGrp)
          annoID <- as.matrix(xsAnnotate@annoID)
          peaklist <- getPeaklist(xsAnnotate)
          rpeaklist <- peaklist[,c("mz","rt","isotopes","adduct","pcgroup")]

	getMi <- function(isoIDidx) {
          molIon <- data.frame(mass=0, retentionTime=0, isotope=0, adduct=0, trace=0)
#          for (i in 1:length(isoID[isoIDidx,1])) {
                 #coord0 <- which(isoID[,1]==isoID[i,1])
                 coord0 <- isoIDidx 
                 molIon0 <- data.frame(matrix(0,length(coord0)+1,5,
                                     dimnames=list(1:(length(coord0)+1),c("mass", "retentionTime", "isotope", "adduct", "trace"))))
		if(polarity=="positive") {
                	 molIon0[1,1] <-  rpeaklist[isoID[coord0[1],1],1]*isoID[coord0[1],4] - (1.007276 * isoID[coord0[1],4])
		}
		if(polarity=="negative") {
                	 molIon0[1,1] <-  rpeaklist[isoID[coord0[1],1],1]*isoID[coord0[1],4] + (1.007276 * isoID[coord0[1],4])
		}
                 molIon0[1,2] <-  rpeaklist[isoID[coord0[1],1],2]
                 molIon0[1,5] <- isoID[coord0[1],1]
		if(polarity=="positive") {
                 	molIon0[-1,1] <- rpeaklist[isoID[coord0,2],1]*isoID[coord0,4] - (1.007276 * isoID[coord0,4])
		}
		if(polarity=="negative") {
                 	molIon0[-1,1] <- rpeaklist[isoID[coord0,2],1]*isoID[coord0,4] + (1.007276 * isoID[coord0,4])
		}
                 molIon0[-1,2] <- rpeaklist[isoID[coord0,2],2]
                 molIon0[-1,3] <- isoID[coord0,3]
                 molIon0[-1,5] <- isoID[coord0,2]
                 molIon <- rbind(molIon, molIon0)
         #}
         molIon <- unique(molIon[-1,])
	return(molIon)
	} 
	molIon <- do.call("rbind", tapply(1:nrow(isoID), isoID[,1], getMi))

         for(i in 1:length(annoGrp[,1])) {
                 molIon0 <- data.frame(matrix(0,1,5, dimnames=list(1,c("mass", "retentionTime", "isotope", "adduct", "trace"))))
                 coord <- which(annoID[,2]==annoGrp[i,1])
                 trace <- annoID[coord,1]
                 molIon0[1,1] <- annoGrp[i,2]
                 molIon0[1,2] <- mean(rpeaklist[trace,2])
                 molIon0[1,4] <-  annoGrp[i,1]
                 molIon0[1,5] <- paste(trace,collapse=";")
                 molIon <- rbind(molIon, molIon0)
         }
	if(allowMiss) {
		anidx <- (unique(unlist(sapply(molIon$trace, function(x) strsplit(x,";")[[1]]))))

		pheno <- xsAnnotate@xcmsSet@phenoData 
		phenon <- rownames(pheno) 
		if(!is.null(toexclude)) {
		phenon <- phenon[-which(pheno[,1] %in% toexclude)] 
		pheno <- pheno[-which(pheno[,1] %in% toexclude),] 
	}
	else {
		pheno <- pheno[,] 
	}

	ntrace <- rownames(peaklist[-as.numeric(anidx),]) 
	totalns <- length(phenon) 

	#peaklist <- getPeaklist(xsAnnotate)

	if(length(unique(as.matrix(pheno)))>1){
	 	ppresence <- apply(peaklist[-as.numeric(anidx),unique(as.matrix(pheno)) ], 1, sum)/totalns 
	}
	else {
	ppresence <- peaklist[-as.numeric(anidx),unique(as.matrix(pheno)) ]/totalns
	}
	cidx <-sub("(^\\d)","X\\1",phenon)
	cidx <- gsub("-", ".", cidx) 
	cidx <- gsub("~", ".", cidx) 
	# experimental
	cidx <- gsub("\\s+", ".", cidx) 
	intmean <- apply(peaklist[-as.numeric(anidx),cidx], 1, function(x) mean(as.matrix(x)))

	smol <- cbind(peaklist[-as.numeric(anidx),c("mz","rt") ][ppresence>minsamp & intmean>minint,],0,0, ntrace[ppresence>minsamp & intmean>minint])
colnames(smol) <- c("mass", "retentionTime", "isotope",  "adduct", "trace") 

	if(polarity=="positive") {
		smol$mass <- smol$mass - 1.007825
	}
	else {
		smol$mass <- smol$mass + 1.007825
	}
	ladduct <- molIon[nrow(molIon), "adduct"]
	smol[,"adduct"] <- ladduct:(nrow(smol)+ladduct-1)
	molIon <- rbind(molIon, smol)
}
	if(sum(duplicated(molIon[,1:2])))	molIon <- molIon[-which(duplicated(molIon[,1:2])),]	
	 molIon$pcgroup <- rpeaklist[as.numeric(sapply(molIon[,"trace"], function(x) strsplit(x, ";")[[1]][1])), "pcgroup"] 

	if(!is.null(xset)) {
		gpidxIso <- as.numeric(molIon[molIon$isotope==1, "trace"])
		gpidx <- as.numeric(molIon[molIon$isotope==0 & molIon$adduct==0, "trace"])
		snList <- list()
		snList$sn <- sapply(gpidxIso, function(x) xset@peaks[xset@groupidx[[x]],"sn"]) 
		snList$sample <- sapply(gpidx, function(x) xset@peaks[xset@groupidx[[x]],"sample"]) 
		snList$sampleIso <- sapply(gpidxIso, function(x) xset@peaks[xset@groupidx[[x]],"sample"]) 
		snList$intb <- sapply(gpidx, function(x) xset@peaks[xset@groupidx[[x]],"intb"]) 
		snList$intbIso <- sapply(gpidxIso, function(x) xset@peaks[xset@groupidx[[x]],"intb"]) 
         	return(list(molIon=molIon,cameraobj=xsAnnotate, snList=snList))
	}
	else {
         return(list(molIon=molIon,cameraobj=xsAnnotate))
	}
}

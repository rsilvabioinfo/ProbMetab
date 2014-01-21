#' get.Mzmatch.annot
#'
#' This function extracts annotation from mzMatch PeakML file,
#' generating a matrix of non-redundant putative molecular ions.
#' 
#' @param filename PeakML file containing ion annotation. 
#' @param onlyBP logical, if TRUE retrieves only PeakML "bp" 
#' relationship, if FALSE also retrieves "potential bp" relationship. 
#' @return A list with a matrix of possible molecular ions with a trace of their annotation, and the used xsAnnotate object. 
#'
#' @export

get.Mzmatch.annot <- function(filename, onlyBP=TRUE) {

nxset <- mzmatch.R::PeakML.xcms.read(filename, annotation=TRUE,version.1=FALSE) 

#cameraobj1 <- xsAnnotate(nxset[[1]]) 
#nxset2 <- fillPeaks(nxset[[1]])  
#cameraobj2 <- xsAnnotate(nxset2) 
#peakl1 <- getPeaklist(cameraobj1) 
#peakl2 <- getPeaklist(cameraobj2) 
#
#idx <- which((peakl1[,"mz"] %in% peakl2[,"mz"]) & (peakl1[,"rt"] %in% peakl2[,"rt"]))
#nxset[[2]] <- nxset[[2]][idx,] 
#
#get.iso.pattern <- function(relationships, onlyBP) {
#selIso <- c(0,0)
#
#for(i in 1:length(relationships)) {
#	if(relationships[i]=="bp"){
#		selIso	<- rbind(selIso, c(i,0))
#	}
#	
#	if(relationships[i]=="bp|C13 isotope #1") {
#		selIso[selIso[,1]==max(which(relationships[1:i]=="bp")),2] <- i
#	}
#	if(!onlyBP) {
#		if(relationships[i]=="potential bp"){
#			selIso	<- rbind(selIso, c(i,0))
#		}
#		if(relationships[i]=="potential bp|C13 isotope #1"){
#			selIso[selIso[,1]==max(which(relationships[1:i]=="potential bp")),2] <- i
#		}
#	}
#}
#return(selIso[-1,])
#}
#
#	isoTable <- get.iso.pattern(nxset[[2]][,3],onlyBP)
#
	rmatrix <- nxset[[2]]
	rmatrix[,1] <- 1:nrow(rmatrix) 
	relationship <- rmatrix[,c(1,3)]
	getBp <- function(idx) {
		relationship <- matrix(relationship[idx,], ncol=2)
		bp <- relationship[which(relationship[,2]==relation),1] 
		bp <- as.numeric(bp) 
		bpc13 <- relationship[which(relationship[,2]==paste(relation, "|C13 isotope #1", sep="")),1]  
		bpc13 <- as.numeric(bpc13) 
		if(length(bp)==length(bpc13)) return(cbind(bp, bpc13)) else return(cbind(bp, c(bpc13,rep(0, length(bp)-length(bpc13)))))
	}

	# warning: this function only holds true if all "C13 isotope #1" corresponds to its firt relationship
	relation <- "potential bp"
	pbp <- do.call("rbind",tapply(1:nrow(rmatrix), rmatrix[,2], getBp)) 

	relation <- "bp"
	bp <- do.call("rbind",tapply(1:nrow(rmatrix), rmatrix[,2], getBp)) 

	isoTable <- rbind(bp, pbp) 
	isoTable <- isoTable[order(isoTable[,1]),] 

	molIon <- rep(0,5); 
	for(i in 1:nrow(isoTable)) { 
		if(isoTable[i,2]==0){ 
			m1 <- matrix(c(nxset[[1]]@groups[isoTable[i,1],c("mzmed", "rtmed")], 0, i, isoTable[i,1]), nrow=1) 
			if(length(m1)!=5) stop("\nSomething wrong\n")
			molIon <- rbind(molIon, m1)
		}
 		else { 
			m1 <- cbind(nxset[[1]]@groups[isoTable[i,],c("mzmed", "rtmed")],c(0,1),0,isoTable[i,])
			if(nrow(m1)!=2) stop("\nSomething wrong\n")
			molIon <- rbind(molIon, m1) 
		}
	} 
 molIon <- molIon[-1,] 
molIonIso <- molIon[molIon[,4]==0,]
sq <- seq(1, nrow(molIonIso), by=2) 

for(j in sq) {
if(abs(diff(molIonIso[j:(j+1),1]))>1.5) {
	srmatrix <- rmatrix[rmatrix[,2]==rmatrix[molIonIso[j,5],2],]
	ridx <- as.numeric(srmatrix[grep("C13 isotope #1", srmatrix[,3]), 1]) 
	val <- nxset[[1]]@groups[ridx, "mzmed"] - molIonIso[j,1]
	nidx <- ridx[which(val>0 & val < 1.5)]
	if(length(nidx)) {
		molIonIso[j+1,1] <- nxset[[1]]@groups[nidx, "mzmed"]
		molIonIso[j+1,5] <- nidx
	}
	else {
		molIonIso[j,4] <- j
		molIonIso[j+1,4] <- -j
	}

}
}
molIonIso <- molIonIso[-which(molIonIso[,4]<0),]
molIon <- rbind(molIonIso, molIonIso[molIonIso[,4]>0,], molIon[molIon[,4]>0,]) 
 
 colnames(molIon) <- c("mass", "retentionTime", "isotope", "adduct", "trace") 
 # for compatibility with get.annot output
 molIon <- as.data.frame(molIon)
 molIon$pcgroup <- 0

#nxset2 <- fillPeaks(nxset[[1]]) 
cameraobj <- xsAnnotate(nxset[[1]])
#cameraobj <- new("xsAnnotate",psSamples=numeric(0), sample=numeric(0), ruleset=data.frame(0), runParallel=list(0), xcmsSet = nxset2)
xset <- nxset[[1]]
	if(!is.null(xset)) {
		gpidxIso <- as.numeric(molIon[molIon$isotope==1, "trace"])
		gpidx <- as.numeric(molIon[molIon$isotope==0 & molIon$adduct==0, "trace"])
		snList <- list()
		snList$sn <- sapply(gpidxIso, function(x) xset@peaks[xset@groupidx[[x]],"sn"]) 
		snList$sample <- sapply(gpidx, function(x) xset@peaks[xset@groupidx[[x]],"sample"]) 
		snList$sampleIso <- sapply(gpidxIso, function(x) xset@peaks[xset@groupidx[[x]],"sample"]) 
		snList$intb <- sapply(gpidx, function(x) xset@peaks[xset@groupidx[[x]],"intb"]) 
		snList$intbIso <- sapply(gpidxIso, function(x) xset@peaks[xset@groupidx[[x]],"intb"]) 
         	return(list(molIon=molIon,cameraobj=cameraobj, snList=snList))
	}
	else {
		return(list(molIon=molIon, cameraobj=cameraobj))
	}

}

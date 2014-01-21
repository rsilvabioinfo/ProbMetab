#' combineMolIon 
#'
#' This function combines ion annotations in different acquisition modes. It operates 
#' in two main modes, combining individual annotations given by get.annot 
#' function, using the retention time and mass/charge windows provided by the user
#' or extracting annotations from a peak table provided by CAMERA's combinexsAnnos
#' function.
#' @param antPOS positive annotation list given by get.annot. 
#' @param antNEG negative annotation list given by get.annot. 
#' @param peaklist given by CAMERA's combinexsAnnos function. If this option
#' is chosen the user has to set the acquisition mode to the same as 
#' in CAMERA's function, and provide the respective object for downstream analysis.
#' @param cameraobj xsAnnotate object for downstream analysis.
#' @param polarity the same CAMERA's function acquisition mode. 
#' @param rtwin retention time window to annotate a peak as present 
#' in both acquisition modes. 
#' @param mzwin mass to charge ratio window to annotate a peak as present 
#' in both acquisition modes. 
#' @return a list with a matrix of possible molecular ions with a 
#' trace of their annotation and the respective xsAnnotate object. 
#'
#' @export

combineMolIon <- function(antPOS, antNEG, peaklist=NULL, cameraobj=NULL, polarity=NULL, rtwin=5, mzwin = 0.05) {
if(!is.null(peaklist)) {
	peakidx <- which(peaklist[,"isotopes"]!="" | peaklist[,"adduct"]!="")
	antPOS3 <- peaklist[peakidx, c("mz", "rt", "isotopes", "adduct")]
	isoidx <- which(antPOS3$isotopes!="") 
	iso <- antPOS3$isotopes[antPOS3$isotopes!=""] 
	 iso <- as.numeric(sapply(iso, function(x) sub("^\\[(\\d+)\\].+", "\\1", x))) 
	charge <- (sapply(antPOS3$isotopes[antPOS3$isotopes!=""][order(iso)], function(x) sub(".+(\\d)\\+|\\-$", "\\1", x))) 
	charge <- suppressWarnings(as.numeric(charge) )
	 charge[is.na(charge)] <- 1 
	niso <- (sapply(antPOS3$isotopes[antPOS3$isotopes!=""][order(iso)], function(x) sub(".+(\\[M.*\\]).+", "\\1", x))) 
	niso <- sub(".+(\\d).+", "\\1", niso) 
	niso <- suppressWarnings(as.numeric(niso)) 
	 niso[is.na(niso)] <- 0 
	preAnt <- cbind(iso[order(iso)], charge, niso, antPOS3[isoidx[order(iso)], c("mz", "rt", "isotopes")])

		  molIon <- data.frame(mass=numeric(nrow(preAnt)), retentionTime=numeric(nrow(preAnt)), isotope=numeric(nrow(preAnt)), adduct=numeric(nrow(preAnt)), trace=numeric(nrow(preAnt)))
			if(polarity=="pos") {
				 molIon$mass <-  preAnt$mz*preAnt$charge - (1.007276 * preAnt$charge)
			}
			if(polarity=="neg") {
				 molIon$mass <-  preAnt$mz*preAnt$charge + (1.007276 * preAnt$charge)
			}
		molIon$retentionTime <- preAnt$rt 
		molIon$isotope <- preAnt$niso
		molIon$adduct <- 0
		molIon$trace <- as.numeric(rownames(antPOS3[isoidx[order(iso)],]))
		molIon$comb <- polarity 

	addidx <- which(antPOS3$adduct!="")  
	v0 <- c(0,0) 
	for(i in 1:length(addidx)) { 
	 	v1 <- suppressWarnings(as.numeric(strsplit(antPOS3$adduct[addidx][i], " ")[[1]]))
		v0 <- rbind(v0, cbind(i, v1[-which(is.na(v1))])) 
	} 

	v0 <- v0[-1,] 
	rnames <- as.numeric(rownames(antPOS3[addidx,])) 
	rnames <-rnames[v0[,1]] 

		  molIon2 <- data.frame(mass=numeric(length(rnames)), retentionTime=numeric(length(rnames)), isotope=numeric(length(rnames)), adduct=numeric(length(rnames)), trace=numeric(length(rnames)))
		
		molIon2$mass <- v0[,2]
		molIon2$retentionTime <- peaklist[rnames, "rt"] 
		molIon2$isotope <- 0  
		molIon2$adduct <- 1:length(rnames) 
		molIon2$trace <- rnames 
		molIon2$comb <- polarity 
		molIon2$comb[which(peaklist[molIon2$trace, ncol(peaklist)]!="")] <- "both" 
		molIon=rbind(molIon, molIon2)
		molIon$pcgroup <- peaklist[as.numeric(sapply(molIon[,"trace"], function(x) strsplit(as.character(x), ";")[[1]][1])), "pcgroup"] 
		
		
	if(sum(duplicated(molIon[,1:2])))	molIon <- molIon[-which(duplicated(molIon[,1:2])),]	
		antComb <- list(molIon=molIon, cameraobj=cameraobj)
		return(antComb)	
}

vidx <- c("",""); 
nvidx <- c("",""); 
pvidx <- c("",""); 

for(i in 1:nrow(antNEG$molIon)) { 
	idx <- which(antPOS$molIon[,1] > (antNEG$molIon[i,1]-mzwin) 
		     & antPOS$molIon[,1] < (antNEG$molIon[i,1]+mzwin) 
		     & antPOS$molIon[,2] > (antNEG$molIon[i,2]-rtwin) 
		     & antPOS$molIon[,2] < (antNEG$molIon[i,2]+rtwin)
		    ) 
	if(length(idx)){
		for(k in 1:length(idx)) {
	# Possible cases
	# same isotopic distribution - discard one
	# adduct convergence - discard one
	# treated as isotope in one, and adduct in another
	#	if the adduct is equal the 12C peak it is discarded
	#	else keep the adduct and the 13C peak 
			if(antNEG$molIon[i,3] != antPOS$molIon[idx[k],3]) {
		 		vidx <- rbind(vidx, c(i, idx[k])) 
				next
			}
			if(antNEG$molIon[i,4] == 0 & antPOS$molIon[idx[k],4]!=0) {
		 		pvidx <- rbind(pvidx, c(i, idx[k])) 
				next
			}
			if(antNEG$molIon[i,4] != 0 & antPOS$molIon[idx[k],4]==0) {
		 		nvidx <- rbind(nvidx, c(i, idx[k])) 
				next
			}
			if(antNEG$molIon[i,4] != 0 & antPOS$molIon[idx[k],4]!=0) {
		 		nvidx <- rbind(nvidx, c(i, idx[k])) 
				next
			}
			if(antNEG$molIon[i,4] == 0 & antPOS$molIon[idx[k],4]==0) {
		 		nvidx <- rbind(nvidx, c(i, idx[k])) 
				next
			}
		 vidx <- rbind(vidx, c(i, idx[k])) 
		}
	}
	
	
}

# for debugging only
 if(!is.null(dim(vidx))) vidx <- vidx[-1,]

  if(!is.null(dim(nvidx))) nvidx <- nvidx[-1,]
  if(!is.null(dim(pvidx))) pvidx <- pvidx[-1,]

antNEG$molIon$ind <- 0
id1 <- which(antNEG$molIon[,3]==1)
antNEG$molIon$ind[id1] <- 1:length(id1)
id1 <- which(antNEG$molIon[,3]==1)-1
antNEG$molIon$ind[id1] <- 1:length(id1)

#antNEG2 <- antNEG$molIon[-as.numeric(vidx[,1]),]
antNEG$molIon$comb <- "neg"
  if(!is.null(dim(pvidx))) antNEG$molIon$comb[as.numeric(pvidx[,1])] <- "both" 

  if(!is.null(dim(nvidx))) {
	 antNEG2 <- as.data.frame(antNEG$molIon[-as.numeric(nvidx[,1]),]) 
  }
  else {
	antNEG2 <- as.data.frame(antNEG$molIon) 
  }

err1 <- setdiff(which(antNEG2[,3]==1), 
		which(antNEG2[,3]==0 & antNEG2[,4]==0)+1

	)

 if(length(err1)) antNEG2 <- as.data.frame(antNEG2[-err1,]) 

err2 <- setdiff(
		which(antNEG2[,3]==0 & antNEG2[,4]==0),
		which(antNEG2[,3]==1) -1
	)

 if(length(err2)) antNEG2 <- as.data.frame(antNEG2[-err2,]) 

sn1 <- unique(antNEG2$ind)
sn1 <- sn1[-which(unique(antNEG2$ind)==0)]
snListNeg <- lapply(antNEG$snList, function(x) x[sn1])

antPOS$molIon$ind <- 0
id1 <- which(antPOS$molIon[,3]==1)
antPOS$molIon$ind[id1] <- 1:length(id1)
id1 <- which(antPOS$molIon[,3]==1)-1
antPOS$molIon$ind[id1] <- 1:length(id1)

antPOS$molIon$comb <- "pos"
  if(!is.null(dim(nvidx))) antPOS$molIon$comb[as.numeric(nvidx[,2])] <- "both" 
  if(!is.null(dim(pvidx))) { 
	antPOS2 <- as.data.frame(antPOS$molIon[-as.numeric(pvidx[,2]),])
  }
  else {
	antPOS2 <- as.data.frame(antPOS$molIon)
  }
err <- setdiff(which(antPOS2[,3]==1), 
		which(antPOS2[,3]==0 & antPOS2[,4]==0)+1
	)
 if(length(err)) antPOS2 <- as.data.frame(antPOS2[-err,]) 

sn1 <- unique(antPOS2$ind)
sn1 <- sn1[-which(unique(antPOS2$ind)==0)]
snListPos <- lapply(antPOS$snList, function(x) x[sn1])

#vpos <- unlist(sapply(vidx[,2], function(x) strsplit(x, ";")[[1]])) 
#antPOS2$comb[as.numeric(vpos)] <- "both" 

antComb <- list(molIon=rbind(antPOS2, antNEG2), antPOS=antPOS, antNEG=antNEG,  
		neg=antNEG$cameraobj, pos=antPOS$cameraobj,
		snListPos=snListPos, snListNeg=snListNeg
	) 
}

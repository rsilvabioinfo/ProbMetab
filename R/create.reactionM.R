#' create.reactionM
#'
#' This function matches a mass vector against a user provided database, 
#' inside an user provided mass tolerance window.
#' @param DB dataframe containing the mandatory fields id, formula 
#' and reactions.
#' @param molIon list of annotations provided by get.annot function. 
#' @param ppm.tol parts per million mass tolerance allowed 
#' in the mass search.
#' @return A matrix of reactions that each compound candidate, inside the
#' mass window, can participate in the metabolism. 
#'
#' @export

create.reactionM <- function(DB, molIon, ppm.tol) {
#	if(is.null(mass.vec)) {
		vmass <- molIon$molIon
		vmass$id <- 1:nrow(vmass)
		vmass <- as.matrix(vmass[vmass[,"isotope"]==0,c("retentionTime","mass", "id")])
#	}
#	else {
#		vmass <- as.matrix(mass.vec)
#	}
	massList <- sapply(DB[,"formula"], function(x) formula2mass(as.character(x)))
	coordList <- lapply(massList, function(x) is.null(x)) 
	massList[unlist(coordList)] <- 0 
	db.mass <- unlist(massList) 
	#ncol <- 4 + max(unlist(lapply(strsplit(as.vector(DB[,5]),";"), length))) 
	reactionM <- matrix("",ncol=6)
	for (i in 1:length(vmass[,2])) {
		mass <- vmass[i,2]
		#logical <- mass < db.mass + mass.tol/2 & mass > db.mass - mass.tol/2
		logical <- abs(((mass-db.mass)/db.mass)*10^6) < ppm.tol		
	if (sum(logical)) {

			#rlist <- strsplit(as.vector(DB[which(logical),5]),";") 
			for (j in 1:length(which(logical))) {
				reactionM0 <- matrix("", nrow=1, ncol=6)
				reactionM0[1,1] <- vmass[i,3] 
				reactionM0[1,2:3] <- vmass[i,1:2]
				reactionM0[1,4] <- db.mass[which(logical)][j]
				reactionM0[1,5] <- as.vector(DB[which(logical),"id"])[j]
				reactionM0[1,6] <- as.vector(DB[which(logical),"reactions"])[j]
#				if(length(rlist[[j]])) {
#					reactionM0[1, 6:(4+length(rlist[[j]]))] <- rlist[[j]] 
#				}
				reactionM <- rbind(reactionM, reactionM0) 
			}
		}
		else {
			reactionM0 <- matrix("", nrow=1, ncol=6)
			reactionM0[1,1] <- vmass[i,3]  
			reactionM0[1,2:3] <- vmass[i,1:2]
			reactionM0[1,5] <- "unknown" 
			reactionM <- rbind(reactionM, reactionM0) 
		}
	}	
	#header <- rep("", ncol(reactionM))
	header <- c("molIonID", "rt", "massObs", "massDB", "id", "reactions")
	colnames(reactionM) <- header
	return(reactionM[-1,])
}

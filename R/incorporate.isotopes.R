#' incorporate.isotopes 
#'
#' Calculates the theoretical pattern of first 13C isotope for each candidate formula.
#'
#' @param ionAnnot annotation list from get.annot function. 
#' @param reactionM compound's reaction matrix.
#' @param comb 1 for acquisition mode combination. 
#' @param polarity acquisition mode polarity. 
#' @param var 1 to use standard mean/sd estimators to carbon number
#' prediction, 2 for median/mad estimators. 
#' @param samp sample indexes, other than blanks, controls and QCs,
#' according to xcms's phenoData.  
#' @param DB data.frame of compound information, with chemical formula. 
#' @return matrix of candidate compound theoretical isotope patterns.
#'
#' @export


incorporate.isotopes <- function(ionAnnot, reactionM, comb=NULL, polarity=NULL,
				 var=2, samp=NULL, DB) {

	#samp <- 12:23 
	#var <- 2



react <- reactionM[reactionM[,"id"]!="unknown",]

isoPatt <- data.frame(id=react[,"id"], formula="", cnun=0, 
		predcnun=0, sd=0, snr=0,
		massObs=as.numeric(react[,"massObs"]), 
		massDB=as.numeric(react[,"massDB"]), 
		isIso=0
	)

isoPatt$formula <- sapply(react[,"id"], function(x) as.character(DB$formula[DB$id==x]))
isoPatt$cnun <- sapply(isoPatt$formula, function(x) sub("C(\\d+)?.+$","\\1",x))
isoPatt$cnun[which(isoPatt$cnun == "")] <- 1 
isoPatt$cnun[!grepl("C", isoPatt$formula)] <- 0 


get.cnun <- function(snList, pos, samp, var) {
	insec <- intersect(snList$sample[[pos]], 
			   snList$sampleIso[[pos]]
		)
	sample <- snList$sample[[pos]]
	sample <- sample[sample %in% insec]
	sample <- sample[sample %in% samp]
	if(!length(sample)) return(rep(NA,3))
	sampleIso <- snList$sampleIso[[pos]]
	sampleIso <- sampleIso[sampleIso %in% insec]
	sampleIso <- sampleIso[sampleIso %in% samp]
	if(!length(sampleIso)) return(rep(NA,3))

	intb <- snList$intb[[pos]][snList$sample[[pos]] %in% sample]
	intbIso <- snList$intbIso[[pos]][snList$sampleIso[[pos]] %in% sampleIso]
	sn <- snList$sn[[pos]][snList$sampleIso[[pos]] %in% sampleIso]

	# max 13C intensity, should change to snr or 12C?
	vid <- tapply(intbIso, sampleIso, which.max) 	
	intv <- tapply(intb, sample, max) 	
	intvIso <- tapply(intbIso, sampleIso, max) 	
	lsn <- tapply(sn, sampleIso, function(x) x ) 	
	snv <- mapply(function(x, y) x[y], lsn, vid) 	
	
	est <- function(x,y) (y*0.9889)/(x*0.0111)
	if(var==1) {
		predcnun <- mean(mapply(est, intv, intvIso)) 
		# consider mad
		sd <- sd(mapply(est, intv, intvIso)) 
		snr <- mean(snv)
		return(c(predcnun, sd, snr))
	}
	else {
		predcnun <- median(mapply(est, intv, intvIso)) 
		# consider mad
		sd <- mad(mapply(est, intv, intvIso)) 
		snr <- median(snv)
		return(c(predcnun, sd, snr))
	}
}

if(is.null(comb)) {
	ionAnnot$molIon$ind <- 0 
	ionAnnot$molIon$ind[which(ionAnnot$molIon[,"isotope"]==1)-1] <- 1:sum( ionAnnot$molIon[,"isotope"]==1)
	ionAnnot$molIon$comb <- polarity 
}
else {
	ind1 <- unique(ionAnnot$molIon$ind[ionAnnot$molIon$isotope==1 & ionAnnot$molIon$comb!="neg"]) 
	ind2 <- unique(ionAnnot$molIon$ind[ionAnnot$molIon$isotope==1 & ionAnnot$molIon$comb=="neg"]) 
}


# match intensities and snrs

for(i in 1:nrow(react)) {

	if(ionAnnot$molIon[as.numeric(react[i,"molIonID"]),"isotope"]==0 &
	   ionAnnot$molIon[as.numeric(react[i,"molIonID"]),"adduct"]==0) {
	if(is.null(comb)) {
		idx1 <- ionAnnot$molIon[as.numeric(react[i,"molIonID"]),c("ind", "comb")]
			pos <- as.numeric(idx1[1])
			if(!is.null(samp)) 1:nrow(ionAnnot$cameraobj@xcmsSet@phenoData)  
			cinf <- get.cnun(ionAnnot$snList, pos, samp, var) 
			isoPatt[i, "predcnun"] <- cinf[1]
			isoPatt[i, "sd"] <- cinf[2]
			isoPatt[i, "snr"] <- cinf[3]
			isoPatt[i, "isIso"] <- 1 
		}
	else {
		idx1 <- ionAnnot$molIon[as.numeric(react[i,"molIonID"]),c("ind", "comb")]
		if(idx1[2]!="neg") {
			pos <- which(ind1==as.numeric(idx1[1]))
			if(!is.null(samp)) 1:nrow(ionAnnot$pos@xcmsSet@phenoData)  
			cinf <- get.cnun(ionAnnot$snListPos, pos, samp, var) 
			isoPatt[i, "predcnun"] <- cinf[1]
			isoPatt[i, "sd"] <- cinf[2]
			isoPatt[i, "snr"] <- cinf[3]
			isoPatt[i, "isIso"] <- 1 
		}
		else {
			pos <- which(ind2==as.numeric(idx1[1]))
			cinf <- get.cnun(ionAnnot$snListNeg, pos, samp, var) 
			isoPatt[i, "predcnun"] <- cinf[1]
			isoPatt[i, "sd"] <- cinf[2]
			isoPatt[i, "snr"] <- cinf[3]
			isoPatt[i, "isIso"] <- 1 
		}
	}
	}
}

return(isoPatt)
}

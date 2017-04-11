#' build.database.kegg
#'
#' Get KEGG compound information needed for
#' exact mass searching and modeling. 
#' @param orgID KEGG's organism id. If NULL the function recovers
#' information from all database compounds. 
#' @return A data.frame with unique id, name, formula
#' and reactions as needed by ProbMetab. 
#'
#' @export


build.database.kegg <- function(orgID=NULL) {
	cpmatrix <- read.delim("http://rest.kegg.jp/list/compound", header=FALSE) 
	rnmatrix <- read.table("http://rest.kegg.jp/link/compound/reaction") 
	cpList <- sapply(1:110, function(n) try(read.delim(paste("http://rest.kegg.jp/find/compound/C", n, "/formula", sep=""), header=FALSE), TRUE)) 
	cpList <- cpList[-grep("Error", cpList)] 
	cpnames <- do.call("rbind", cpList)
	if(!is.null(orgID)) {
		#orgSpecific <- read.table(paste("http://rest.kegg.jp/link/", orgID, "/compound", sep=""))
		orgSpecific <- read.table(paste0("http://rest.kegg.jp/link/", orgID, "/pathway"))
		allpath <- read.table("http://rest.kegg.jp/link/pathway/cpd")
		orgSpecific <- allpath[sub(".+(\\d{5}$)", "\\1",  allpath[,2]) %in% sub(".+(\\d{5}$)", "\\1",  orgSpecific[,1]),]
		cpnames <- cpnames[which(cpnames[,1]%in%orgSpecific[,1]),] 
		missing <- setdiff(unique(orgSpecific[,1]), cpnames[,1]) 
	}
	else {
		missing <- setdiff(cpmatrix[,1], cpnames[,1])
	}
	if(length(missing)) {
		missformula <- sapply(missing, get.formula.kegg) 
		if(sum(which(unlist(lapply(missformula, length))==0))) {
			missing <- missing[-which(unlist(lapply(missformula, length))==0)] 
			missformula <- missformula[-which(unlist(lapply(missformula, length))==0)]
			cpformulas <- rbind(as.matrix(cpnames), cbind(missing, unlist(missformula)))
		}
	}
	else {
		cpformulas <- cpnames
	}
	rn <- sapply(cpmatrix[,1], function(x) paste(rnmatrix[which(rnmatrix[,2]==as.character(x)),1], collapse=";"))
	names <- sapply(cpmatrix[,2], function(x) strsplit(as.character(x), ";")[[1]][1])

	cpformulatotal <- data.frame(id=cpformulas[,1], formula=cpformulas[,2]) 
	cptotal <- data.frame(id=cpmatrix[,1], name=names, reactions=rn)
	cptotal <- merge(cptotal, cpformulatotal, by.x="id", by.y="id") 
	return(cptotal)
}

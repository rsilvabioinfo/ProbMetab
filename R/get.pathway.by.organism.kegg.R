#' get.pathway.by.organism.kegg
#'
#' Get KEGG pathway codes from its API. 
#' @param organismId KEGG's organism id. 
#' @return A vector of KEGG pathway codes. 
#'
#' @export



get.pathway.by.organism.kegg <- function(organismId) {

	path <- try(read.table(paste("http://rest.kegg.jp/link/pathway/",organismId, sep="")), TRUE)
	if(!sum(grepl( "Error", path)))	unique(sub("path:\\D+","",path[,2]))
}

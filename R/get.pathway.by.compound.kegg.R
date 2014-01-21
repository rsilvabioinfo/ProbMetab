#' get.pathway.by.compound.kegg
#'
#' Takes KEGG id code and retrieves the pathways
#' that the compound is possible involved. 
#' @param keggId KEGG's compound id. 
#' @return A string with the pathways collapsed by ";" character. 
#'
#' @export



get.pathway.by.compound.kegg <- function(keggId) {

	path <- try(read.table(paste("http://rest.kegg.jp/link/pathway/",keggId, sep="")), TRUE)
	if(!grepl( "Error", path))	paste(unique(sub("path:\\D+","",path[,2])), collapse=";")
}

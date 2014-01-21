#' get.reactions.by.compound.biocyc
#'
#' Get biocyc reaction codes from its API. 
#' @param cpdID code. 
#' @return character vector of reaction concatenated by ";" character.
#'
#' @export

get.reactions.by.compound.biocyc <- function(cpdID) {
# require(XML)
doc <- readLines(paste("http://websvc.biocyc.org/apixml?fn=reactions-of-compound&id=", cpdID, "&detail=[none|low|full]", sep=""))
doc <- xmlTreeParse(doc)
lreac <- xmlSApply(doc, xmlAttrs) 
vreac <- unlist(lapply(lreac, function(x) x["ID"])) 
vcpdReac <- paste(vreac, collapse=";")
return(vcpdReac)
}


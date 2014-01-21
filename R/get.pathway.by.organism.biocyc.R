#' get.pathway.by.organism.biocyc
#'
#' Get biocyc pathway codes from its API. 
#' @param orgID code. 
#' @return character vector of pathway codes. 
#'
#' @export

get.pathway.by.organism.biocyc <- function(orgID) {
# require(XML)
doc <- readLines(paste("http://websvc.biocyc.org/xmlquery?[x:x<-", orgID,"^^pathways]", sep=""))
doc <- xmlTreeParse(doc)
lpth <- xmlSApply(doc, xmlAttrs) 
vpth <- unlist(lapply(lpth, function(x) x["ID"])) 
return(vpth)
}


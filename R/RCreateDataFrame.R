#' RCreateDataFrame 
#'
#' Create a data frame with compound information from
#' mzMatch open source project format.
#' @param file xml input file. 
#' @return a data frame with compound's id, name and formula. 
#'
#' @export

RCreateDataFrame <-
function(file) {
#    require(XML)
    doc <- xmlParse(file) 
    id <- unlist(xpathApply(doc, "/compounds/compound/id",xmlValue)) 
    name <- unlist(xpathApply(doc, "/compounds/compound/name",xmlValue)) 
    formula <- unlist(xpathApply(doc, "/compounds/compound/formula",xmlValue)) 
    formula <- gsub(".+;?\\[|\\].+","",formula)
    data.frame(id,  name, formula)
}

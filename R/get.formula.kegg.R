#' get.formula.kegg
#'
#' Get KEGG compound formulas from its API. 
#' @param x KEGG compound code. 
#' @return character formula. 
#'
#' @export

get.formula.kegg <-
function(x){
        inf <- readLines(paste("http://rest.kegg.jp/get/", x, sep=""))
        coord <- grep("FORMULA", inf) 
        name <- gsub("^FORMULA\\s+|;","",inf[coord])  
        return(name)
}

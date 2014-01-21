#' get.name 
#'
#' Get KEGG object names from its API. 
#' @param x KEGG code. 
#' @return character name. 
#'
#' @export

get.name <-
function(x){
        inf <- readLines(paste("http://rest.kegg.jp/get/", x, sep=""))
        coord <- grep("NAME", inf) 
        name <- gsub("^NAME\\s+|;","",inf[coord])  
        return(name)
}

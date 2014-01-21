#' get.ncbi.id 
#'
#' Transtates KEGG's enzyme code to NCBI gene id. 
#' @param x KEGG's enzime code. 
#' @return a character code of NCBI gene entry. 
#'
#' @export

get.ncbi.id <-
function(x) {
    if(length(readLines(paste("http://rest.kegg.jp/link/genes/", x, sep="")))>1 ){
             x2 <- as.character(read.delim(paste("http://rest.kegg.jp/link/genes/", x, sep=""), header=FALSE)[1,2]) 
             as.character(read.delim(paste("http://rest.kegg.jp/conv/ncbi-gi/",x2,sep=""), header=FALSE)[1,2]) 
     }
}

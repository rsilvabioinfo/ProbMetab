#' get.compounds.by.pathway.biocyc
#'
#' Get biocyc compound codes from its API. 
#' @param pathID code. 
#' @return character matrix of compound codes, names and formulas. 
#'
#' @export

get.compounds.by.pathway.biocyc <- function(pathID) {
# require(XML)
doc <- readLines(paste("http://websvc.biocyc.org/apixml?fn=compounds-of-pathway&id=", pathID, "&detail=[none|low|full]", sep=""))
doc <- xmlTreeParse(doc)
lcpd <- xmlSApply(doc, xmlAttrs) 
vcpd <- unlist(lapply(lcpd, function(x) x["ID"])) 

cloc <- xmlSApply(doc, xmlName)=="Compound" 
seq <- 1:length(xmlRoot(doc)) 
biocyc <- rep("", 3)
for(i in seq[cloc]){
 vcpd <- lcpd[[i]]["ID"]  
 vform <- try(xmlSApply(xmlRoot(doc)[[i]][["cml"]][["molecule"]], xmlAttrs)$formula,TRUE) 
 vname <- try(xmlSApply(xmlRoot(doc)[[i]][["cml"]], xmlAttrs)["title","molecule"], TRUE) 
 biocyc <- rbind(biocyc, c(vcpd, vname, vform))
}
biocyc[,2] <- gsub("\\s+", "" ,biocyc[,2])
#biocyc <- cbind(vcpd, vname, vform)
#suppressWarnings(biocyc <- cbind(vcpd, vname, vform)) 
biocyc <- matrix(biocyc[-1,], ncol=3)
colnames(biocyc) <- c("id", "name", "formula")
row.names(biocyc) <- rep("", nrow(biocyc)) 
return(biocyc)
}


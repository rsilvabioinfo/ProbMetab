#' get.kgml.positions.kegg
#'
#' Gets KEGG's pathway code, download the kgml file, and
#' retrieves the pathway layout. 
#' @param path  KEGG's pathway code.
#' @return a list with an adjacency matrix for pathway nodes,
#' and a matrix with node positions. 
#'
#' @export


get.kgml.positions.kegg <- function(path) {

#require(XML) 
link <- paste("http://www.kegg.jp/kegg-bin/download?entry=", path, "&format=kgml", sep="")

doc <- xmlParse(link,ignoreBlanks = TRUE)  
entry <- getNodeSet(doc, "//entry") 
id <- lapply(entry, xmlGetAttr, "id") 
name <- lapply(entry, xmlGetAttr, "name") 
name <- sub("^\\S*(\\w{1}\\d{5}).*", "\\1", name) 
type <- lapply(entry, xmlGetAttr, "type") 

#id <- xpathApply(doc, "//entry[@id]", xmlGetAttr, "id") 
#name <- xpathApply(doc, "//entry[@name]", xmlGetAttr, "name") 
#type <- xpathApply(doc, "//entry[@type]", xmlGetAttr, "type") 
#
graphics <- getNodeSet(doc, "//graphics")   
#grname <- xpathApply(doc, "//graphics[@name]", xmlGetAttr, "name") 
#grx <- xpathApply(doc, "//graphics[@x]", xmlGetAttr, "x") 
#gry <- xpathApply(doc, "//graphics[@y]", xmlGetAttr, "y") 
#
grx <- lapply(graphics, xmlGetAttr, "x") 
gry <- lapply(graphics, xmlGetAttr, "y")
grx[unlist(lapply(grx, is.null))] <- "" 
gry[unlist(lapply(gry, is.null))] <- ""
grname <- lapply(graphics, xmlGetAttr, "name")
grname <- sub("^\\S*(\\w{1}\\d{5}).*", "\\1", grname) 

cp <- data.frame(id=unlist(id), name=unlist(name), type=unlist(type))
gr <- data.frame(grname=unlist(grname), grx=unlist(grx), gry=unlist(gry))

cp <- merge(cp, gr, by.x="name", by.y="grname")
hcp <- as.matrix(cp[which(cp[,"type"]=="compound"),])

Reacts <- getNodeSet(doc, "//reaction")  

sepReac <- function(rn) {
rnId <- xmlGetAttr(rn, "id") 
reId <- lapply(xmlChildren(rn), xmlGetAttr, "id")
reName <- lapply(xmlChildren(rn), xmlGetAttr, "name")
cbind(rnId, reId, reName) 
}
rnTab <- do.call("rbind", lapply(Reacts, sepReac)) 
rnTab[,3] <- sub("^\\S*(\\w{1}\\d{5}).*", "\\1",rnTab[,3]) 

adj <- matrix(0, nrow=nrow(hcp), ncol=nrow(hcp))
colnames(adj) <- rownames(adj) <- as.matrix(hcp[,"name"]) 

for(i in 1:nrow(hcp)) {
r1 <- rnTab[which(rnTab[,2]%in%hcp[i,"id"]),1] 
for(j in 1:length(r1)) {
idx <- which(hcp[,"name"] %in% rnTab[which(rnTab[,1]%in%r1[j]),3])
adj[idx, idx] <- 1
}
}
 diag(adj) <- 0 

hcp[,"name"] <- sapply(hcp[,"name"], function(x) strsplit(as.character(x), " ")[[1]][1])

dup <- which(duplicated(hcp[,"name"]))
hcp[dup, "name"] <- paste(hcp[dup, "name"], 1:length(dup), sep="")


colnames(adj) <- rownames(adj) <- as.matrix(hcp[,"name"]) 

return(list(adj=adj, posMatrix=hcp))

}

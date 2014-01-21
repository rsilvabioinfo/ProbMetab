#' create.pathway.node.attributes
#'
#' This function writes a standard Cytoscape node attribute list (http://www.cytoscape.org/) file.
#' It takes the compound codes and retrieve all the known pathways,
#' where the compound is known to be present. It only works for KEGG,
#' but, a specification of database will be available soon. 
#' @param classTable classification table created by export.class.table function.  
#' @param graph graphNEL object, with node indexes corresponding to mass indexes in classTable. 
#' @param DB database with compound names associated to unique ids, used by export.class.table function.  
#' @param filename1 filename to attribute pathway list file. 
#' @param filename2 optional filename to attribute pathway list discriminating compound/pathway associations. 
#' @param organismId KEGG organism id (http://www.kegg.jp/kegg/catalog/org_list.html)
#' to filter possibibly pathwyas for known pathways for that organism. 
#' @return writes a standard Cytoscape attribute list to current working directory. Also creates
#' a matrix containing putative compound counting to each pathway. 
#'
#' @export

create.pathway.node.attributes <- function(classTable, graph, DB, filename1, filename2=NULL, organismId=NULL) { 

	for(i in 1:nrow(classTable)) if(classTable[i,1]=="") classTable[i,c(1,4,6,7)] <- classTable[i-1,c(1,4,6,7)]
	msel <- as.matrix(classTable[,1:7])
	msel <- cbind(msel[,6], msel[,-6])
	colnames(msel)[1] <- "Id"
	msel[,1] <- as.numeric(sub("^\\s+","", msel[,1]))

	# a function to generate node and edge attribute lists
	subSel <- sapply(nodes(graph), function(x) msel[which(msel[,1]==x),3]) 

	#library(doMC)
	#registerDoMC()
	#getPList <- function(i) unlist(sapply(subSel[[i]], get.pathway.by.compound.kegg)) 
	#pts <- foreach(i=1:length(subSel)) %dopar% getPList(i)
	#

	cpdPth <- read.table("http://rest.kegg.jp/link/pathway/compound") 

	getPath <- function(cplist) {
		pt1 <- sapply(cplist, function(x) paste(unique(sub("^\\D+", "", cpdPth[grep(x, cpdPth[,1]),2])), collapse=";"))
		pt1 <-  pt1[pt1!=""] 
		return(pt1)
	}

	pts <- lapply(subSel, getPath)
	pts2 <- lapply(pts, function(x) paste(unique(unlist(lapply(x, function(x) strsplit(x, ";")))), collapse=";")) 
	pts2[(unlist(lapply(pts2, is.null)))] <- ""
	pts2 <- unlist(pts2) 
	t1 <-table(strsplit(paste(pts2, collapse=";"),";")[[1]])
	t1 <- t1[-1]
	t1 <- t1[t1>4]

	pnames <- sapply(names(t1), function(x) get.name(paste("path:map", x, sep="")))

	cpdInfo <- cbind(names(t1), pnames, t1)
	toExclude <- which(pts2=="")


	if(!is.null(organismId)) {
		spec <- get.pathway.by.organism.kegg(organismId)
		cpdInfo <- cpdInfo[cpdInfo[,1] %in% spec,] 
		for(i in 1:length(pts)) {
			l1 <- sapply(pts[[i]], function(x) strsplit(x, ";")[[1]] %in% spec)
			tmp <- mapply(function(x, y) paste(strsplit(x, ";")[[1]][y], collapse=";"), pts[[i]], l1)
			pts[[i]] <- tmp[tmp!=""] 
		}
	}

	cpdInfo <- rbind(c("Id", "Pathway Name", "N. Compounds"), cpdInfo)
	#ex <- c(54, 55, 73:85, 87:109)
	#cpdInfo <- cpdInfo[-ex,]

	getIndAttr <- function(list) {
		n <- as.character(sapply(names(list), function(x) DB$name[which(DB$id==x)]))
		paste("(", paste(paste(n, list), collapse="::"), ")", sep="") 
	}

	write.table("Pathways (class=java.lang.String)", filename1, row.names=FALSE, quote=FALSE, col.names=FALSE)
	write.table(paste(nodes(graph)[-toExclude], " = (", gsub(";","::", pts2[-toExclude]),")", sep=""), filename1, row.names=FALSE, col.names=FALSE, quote=FALSE, append=TRUE)
	if(!is.null(filename2)) {
		empty <- which(unlist(lapply(pts, length))==0) 
		write.table("PathwaysDetails (class=java.lang.String)", filename2, row.names=FALSE, quote=FALSE, col.names=FALSE)
		write.table(unlist(paste(nodes(graph)[-empty], "=", lapply(pts, getIndAttr)[-empty])), filename2, row.names=FALSE, col.names=FALSE, quote=FALSE, append=TRUE)
	}
	return(cpdInfo)
}

#' create.reaction.edge.attributes
#'
#' This function writes a standard Cytoscape edge attribute list (http://www.cytoscape.org/) file.
#' It takes the the possible mass connections, codified w matrix, retrieve the reactions
#' where the compound is known to be present and associate then to a mass edge. It only works for KEGG,
#' but, a specification of database will be available soon. 
#' @param classTable classification table created by export.class.table function.  
#' @param graph graphNEL object, with node indexes corresponding to mass indexes in classTable. 
#' @param w matrix of compound connections.
#' @param reactionM data.frame with compound annotation information. 
#' @param DB database with compound names associated to unique ids, used by export.class.table function.  
#' @param filename reaction attribute list file. 
#' @return Writes a standard Cytoscape attribute list to the current working directory. 
#'
#'
#' @export

create.reaction.edge.attributes <- function(classTable, graph, w, reactionM, DB, filename) { 

	mw <- which(w==1, arr.ind=TRUE)
	enames <- edgeNames(graph) 
	rclass <- classTable[classTable[,2]!="unknown",]
	for(i in 1:nrow(rclass)) if(rclass[i,6]=="") rclass[i,6] <- rclass[i-1,6] 
	rclass[,6] <- as.numeric(gsub("\\s+","",rclass[,6]))
	rreact <- reactionM[,-1][reactionM[,-1][,4]!="unknown",] 

	edgeAtt <- c("", "")
	for(j in 1:length(enames)) {
		edgetmp <- ""
		ns <- as.numeric(strsplit(enames[j], "~")[[1]])
		v1 <- mw[,1] %in% which(rclass[,6]==ns[1])
		v2 <- mw[,2] %in% which(rclass[,6]==ns[2])
		cidx <- matrix(mw[which(v1 & v2),], ncol=2) 
		for(i in 1:nrow(cidx)) {
			 vr <- sapply(cidx[i,], function(y) rreact[y,5])
			 r <-  strsplit(vr[1], ";")[[1]][which(strsplit(vr[1], ";")[[1]] %in% strsplit(vr[2], ";")[[1]])] 
			 vc <- sapply(cidx[i,], function(y) rreact[y,4])
			 #vc <- sapply(vc, get.name) 
			 vc <- sapply(vc, function(x) as.character(DB$name[which(DB$id==x)]))
			 edgetmp <- rbind(edgetmp,
					c(paste(vc, collapse=paste(" ", r, " ", sep="")))
				    )
		}
		edgeAtt <- rbind(edgeAtt, c(paste(ns[1], "(unspecified)", ns[2]),
					   paste("(",paste(edgetmp[-1], collapse="::"), ")", sep="")
			   		)
			)
		}
	edgeAtt <- edgeAtt[-1,] 
	write.table("Reactions (class=java.lang.String)", filename, row.names=FALSE, col.names=FALSE, quote=FALSE) 
	write.table(paste(edgeAtt[,1], "=", edgeAtt[,2]), filename, row.names=FALSE, col.names=FALSE, quote=FALSE, append=TRUE) 
}

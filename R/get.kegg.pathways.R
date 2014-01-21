#' get.kegg.pathways 
#'
#' Get pathway map links to all kegg pathways that have more than
#' numCpds as nodes.
#'
#' @param cpds KEGG compound ids, including "cpd:". 
#' @param numCpds minimum number of compounds as nodes of the pathway.
#' 
#' @return list of the missing compounds of kegg pathways, and a table with
#' map links and statistics. 
#'
#' @export



get.kegg.pathways <- function(cpds, numCpds) {

		# find pathways with KEGG rest
		cpds <- sub("\\.",":", cpds) 
		keggPath <- "http://rest.kegg.jp/link/pathway/"
		#cpds <- as.character(colnames(W))
		paths <- c("","") 
		for(i in 1:length(cpds)) {
		paths2 <- try(read.table(paste(keggPath,cpds[i] , sep=""), header=FALSE),TRUE)
			
			if(sum(is(paths2)!="try-error")) paths <- rbind(paths, paths2)
		}
		paths <- paths[-1,]
		missing <- setdiff(cpds, unique(paths[,1]))	
#		missing <- setdiff(colnames(W), unique(paths[,1]))
#		paths <- c("","")
#		while(sum(is(paths2)!="try-error")) {
#			paths2 <- try(read.table(paste(keggPath, paste(missing, collapse="+"), sep=""), header=FALSE), TRUE) 
#			if(sum(is(paths2)!="try-error")) paths <- rbind(paths, paths2)
#			missing <- setdiff(colnames(W), unique(paths[,1]))
#		}

		# produce a unique graph with all compound connections (really? have you tested with reaction bank)
		
		paths[,2] <- gsub("\\D","",paths[,2]) 
		paths <- unique(paths) 
		paths[,2] <- paste("map",paths[,2],sep="")  
		allPath <- unique(paths[,2]) 
		keggMapLink <- "http://www.kegg.jp/kegg-bin/show_pathway?" 	
	mpLinks <- matrix(c("","","","",""), ncol=5)
	for(i in 1:length(allPath)) {		
		pcoord <- grep(allPath[i], paths[,2]) 
		if(length(pcoord) > numCpds) {
		name <- readLines(paste("http://rest.kegg.jp/get/",allPath[i],sep="")) 
		name <- sub("NAME\\s+","",name[grep("NAME", name)])
		mpLinks <- rbind(mpLinks, c(name,allPath[i],length(pcoord), paste(c(paste(keggMapLink, allPath[i], "/default%3dred", sep=""), as.character(paths[pcoord,1])), collapse="/"),paste(as.character(paths[pcoord,1]), collapse=";")))
		}	
	}
	mpLinks <- matrix(mpLinks[-1,], ncol=5)

 	if(length(mpLinks)==0) { return(cat("\nNo free lunch for you mister\n")) } 
	colnames(mpLinks) <- c("Name", "Pathway","Num Compounds","Link","Compounds")
	rownames(mpLinks) <- 1:nrow(mpLinks)
	return(list(mpLinks=mpLinks, missing=missing))
}

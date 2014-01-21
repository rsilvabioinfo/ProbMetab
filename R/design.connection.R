#' design.connection 
#'
#' Design a connection matrix w from unique reaction identifiers
#' from experimentally known reactions between candidate compound.
#'
#' @param reactionM a data.frame with compounds and its connections.
#' @return A binary matrix w of connections between candidate compound. 
#'
#' @export

design.connection <- function(reactionM) {
	ReactMatrix <-reactionM[reactionM[,"id"]!="unknown",]
	w <- matrix(0, nrow(ReactMatrix), nrow(ReactMatrix))

	count <- 0

	connect <- function(y) {
		reacs <- strsplit(as.character(y[6]), ";")[[1]]
		coord <- sapply(reacs, function(x) grep(x, ReactMatrix[,6]))
		coord <- unique(as.vector(unlist(coord)))
		# The candidates of the same mass can't have candidate compounds connected, neither
		# should same compounds assigned to different masses
		rule <- y[5]!=ReactMatrix[coord,5] & y[3]!=ReactMatrix[coord,3]
		coord <- coord[rule]
		count <<- count + 1
		w[count, coord] <<- 1  
	}

	apply(ReactMatrix, 1, connect)
	return(w)
}

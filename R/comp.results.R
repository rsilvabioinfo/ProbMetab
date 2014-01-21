#' comp.results 
#'
#' Compare two classification tables, given by export.class.table,
#' and reports the difference between two different models.  
#'
#' @param reactionM matrix with reactions of each candidate compound. 
#' @param w matrix of compound connections. 
#' @param ansLik a list of mass to compound assignment, 
#' based only on a likelihood.
#' @param ansConn list of mass to compound assignment,
#' with compound's connections contribution.   
#' @return a list with compound classification classes, table index of classes and
#' a matrix of intensities of selected compounds. 
#'
#' @export

comp.results <- function(reactionM, w, ansLik, ansConn) {

	# find a better way to work with redundancies,
	# the w1,2,3 variables suffer from this...
	ansLik <- ansLik$classTable 
	ansConn <- ansConn$classTable 
 	
 	ansLik <- ansLik[ansLik[,2]!="unknown",] 
 	ansConn <- ansConn[ansConn[,2]!="unknown",] 
	kmatrix <- reactionM[reactionM[,4]!="unknown",1:4]
	unq <- unique(kmatrix[,1:2]) 
	connM <- cbind(reactionM[reactionM[,4]!="unknown",4], apply(w, 1, sum))

	# what is the criteria to count as a change for probability?
	# first: likelihood alone does not distinguish candidates, and added
	# information does.
	# second: the criteria for added information is evaluated, as bigger than
	# the remaining candidates... 

	# walk through ans
	v1 <- c() # only one candidate
	v2 <- c() # likelihood enough 
	v3 <- c() # reaction determines 
	v4 <- c() # probability change not determined 
	for (i in 1:nrow(unq)) {
		coord <- which(kmatrix[,1]==unq[i,1] & kmatrix[,2]==unq[i,2])
		if (length(coord)>1) {
			if(ansLik[coord,3][1]=="undef" | ansLik[coord,3][2]=="undef")  { 
			v4 <- append(v4, coord[1])
			next
			}
			if(as.numeric(ansLik[coord,3][1])>as.numeric(ansLik[coord,3][2]) &
			   ansLik[coord,2][1]==ansConn[coord,2][1]) {
				v2 <- append(v2, coord[1])
				next
			}
			else {  
				if(as.numeric(ansConn[coord,3][1])>as.numeric(ansConn[coord,3][2])) {
					w1 <- which(connM[coord,1]==ansConn[coord,2][1])
					w2 <- which(connM[coord,1]==ansConn[coord,2][2])
					w3 <- which(connM[coord,1]==ansLik[coord,2][1])
					if(!length(w1)){
						v4 <- append(v4, coord[1])
						next
					}
					if((as.numeric(connM[coord,2][w1])>as.numeric(connM[coord,2][w2])) | (as.numeric(connM[coord,2][w1])>as.numeric(connM[coord,2][w3]))){ 
						v3 <- append(v3, coord[1])
					}
					else {
						v4 <- append(v4, coord[1])
					}
				}
				else {
			
					v4 <- append(v4, coord[1])
				}
			}			
		}
		else {

			v1 <- append(v1, coord[1])
			next
		}	 
	}
	coList <- list(v1=v1,v2=v2,v3=v3,v4=v4)
	idx <- as.vector(unlist(coList[1:3]))
	W <- w 
	colnames(W) <- rownames(W) <- connM[,1] 
	cnames <- which(colnames(W) %in% ansConn[idx,2])
	W <- W[cnames, cnames]
	dup <- which(duplicated(colnames(W)))
	W <- W[-dup, -dup]
	intMatrix <- ansConn[idx,c(2,8:ncol(ansConn))]
	# I have to change classTable from make.matrix to have
	# colnames
	
	return(list(W=W, idxClass=idx, classList=coList, connM=connM, intMatrix=intMatrix))
}

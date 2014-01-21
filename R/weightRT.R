#' weightRT 
#'
#' Builds a c (number of compounds) by m (number of masses) matrix of weights 
#' 
#' @param rtObj is a list of fitted values predicted by rt.predict. 
#' @param reactionM is the KEGG search list. 
#' @param userCuttoff is the user set expectation of error acceptance.
#' @param rtWeigth is penalty to errors that fall bellow the threshold given b user. 
#' @param plot logical, wheter or not do plot the output of the function. 
#' @return A matrix wrt of likelihood. 
#'
#' @export

weightRT <-
function(rtObj, reactionM, userCuttoff=.95, rtWeigth=0.1, plot=FALSE) {       

	voidTime <- rtObj$voidTime
	predError <- (((10^rtObj$ans$fitted.values)*voidTime+voidTime)-(rtObj$testSet[,2]*voidTime+voidTime))/(rtObj$testSet[,2]*voidTime+voidTime)

	qnt <- quantile(abs(predError), userCuttoff)                      
    #qnt <- quantile(abs(rtObj$ans$predError)/rtObj$testSet$data, userCuttoff)                      
    #realDataqt <- abs(rtObj$pred-rtObj$predSet$data)/rtObj$predSet$data
	realDataqt <- abs( (((10^rtObj$pred)*voidTime+voidTime)-(rtObj$predSet[,2]*voidTime+voidTime))/(rtObj$predSet[,2]*voidTime+voidTime) )
    if(plot) {
    abvCutoff <- dexp(realDataqt[realDataqt<=qnt],qnt)
    blwCutoff <- rtWeigth*dexp(realDataqt[realDataqt>qnt],qnt) 
    plot(realDataqt, c(abvCutoff, blwCutoff), type="n",xlab="Retention Time Error", ylab="Weigth associated to error")  
    points(realDataqt[realDataqt<=qnt], abvCutoff, col=3)                                                                                                      
    points(realDataqt[realDataqt>=qnt], blwCutoff, col=2)                                                                                                      
    abline(v=qnt, col=2)
    }

    wv <- sapply(realDataqt, function(x) if(x<=qnt) dexp(x,qnt) else rtWeigth*dexp(x,qnt))
	reactMatrix <- reactionM[reactionM[,5]!="unknown",]
	o <- (sapply(names(table(reactMatrix[,3])), function(x)which(reactMatrix[,3]==x)[1]))
	t <- table(reactMatrix[,3])[order(o)]
	fac <- rep(1:nrow(t), as.vector(t))
	for(i in 1:length(t)) if(sum(rtObj$pred[fac==i]==rtObj$ans$coef[1])) wv[fac==i] <- 1 
    wm <- matrix(0,nrow=nrow(reactMatrix), ncol=length(t))
    for(i in 1:length(t)) wm[fac==i, i] <- wv[fac==i] 

#    keggObj <- keggObj$m1                                                                                                    
#    #count <- unique(keggObj[,1])
#    count <- unique(keggObj[,c(1,2)]) 
#    count <- as.matrix(count[,1]) 
#    wm <- matrix(0,nrow=nrow(keggObj), ncol=length(count))
# 
#        for(i in 1:length(count)) {
#              coord <- which(keggObj[,1]==count[i])   
#            for (j in 1:length(coord)) {
#                if(length(c <- grep(keggObj[coord[j],4], rtObj$predSet$id))){ # there where molecular descriptors for   
#                    for (k in 1:length(c)) {                                  # with the same retention time
#                        if (rtObj$predSet$data[c[k]]==keggObj[coord[j],1]) { 
#                            wm[coord[j],i] <- wv[c[k]]    
#                         }
#                     }
#                }
#                else if (length(c2 <- grep(count[i], rtObj$predSet$data))){ # Descriptor only for some candidates 
#                    wm[coord[j],i] <- runif(1, min(wv[c2]), max(wv[c2]))    
#                }
#                else { # No descriptors 
#                    wm[coord[j],i] <- 1
#                }
#            }
#        }
    return(list(wm=wm, rtProb=wv))
}

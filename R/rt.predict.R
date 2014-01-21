#' rt.predict 
#'
#' Predict retention time from a given training set,  
#' and a matrix of compound descriptors.
#'
#' @param testSet is a training set provided by user. 
#' @param predSet is a set of candidate compound to have their times predicted. 
#' @param descData is the user provided list of matrices of compound descriptors. 
#' @param descList is the optional user provided list of descriptors name calculated by rcdk package. 
#' @param columnPar is the optional user provided list of column parameters. 
#' @param voidTime column "dead" time. 
#' @return A list of prediction parameters.
#'
#' @export


rt.predict <-
function (testSet, predSet, descData=NULL, descList=NULL, columnPar=NULL, voidTime=NULL) {
#    require(rcdk)
    require(mgcv)
    require(leaps)
    require(bootstrap)
#    require(RCurl)
#    require(foreach)
#    require(doMC) 

    if (!is.null(columnPar)) {
        voidTime <- function(v) {
            # c(d, L, F)
            # d and L in cm
            # F in ml/min 
            # t0 in min
                t0 <- (1/7)*(v[1]^2)*pi*v[2]/v[3] 
                return(t0)
        }
        rf <-function(RT, voidtime) { (RT-voidtime)/voidtime }
        testSet$data <- rf(testSet$data, voidTime(columnPar))
        predSet$data <- rf(predSet$data, voidTime(columnPar))
    }
	
    if (!is.null(voidTime)) {
        rf <-function(RT, voidtime) { (RT-voidtime)/voidtime }
        testSet$data <- rf(testSet$data, voidTime)
        predSet$data <- rf(predSet$data, voidTime)
	}
    else {
        testSet$data <- testSet$data
        predSet$data <- predSet$data
    }
# 10^pred[1]*voidtime -voidtime

# should use
#
# library(SSOAP)
# library(rpubchem)
# kegg <- processWSDL("http://soap.genome.jp/KEGG.wsdl")
# iface <- genSOAPClientInterface(def = kegg)
# iface@functions$get_linkdb_between_databases("cpd:C04443", "PubChem", 1, 10)
# get.sid()

#get.smile1 <- function(keggList) {
#    db <- read.delim("data/ac2021823_si_003.csv") 
#    mywhich <- function(x, db) which(db==x) 
#    l <- sapply(keggList, mywhich, db[,"KEGGid."]) 
#    mylist <- function(l) l[[1]][1] 
#    l <- lapply(l, mylist) 
#    v1 <- which(unlist(lapply(l, length))!=0)
#    sml <- as.vector(db[unlist(l[v1]),20]) 
#    if (sum(sml=="")>0) {
#        v2 <- which(sml=="") 
#        for ( i in 1:length(v2) ) {
#           if( length(get.smile2(keggList[v2[i]])) ) sml[v2[i]] <- get.smile2(keggList[v2[i]])
#           else sml <- sml[-v2[i]] 
#        }
#    }
#    # exclude element when there is no smile
#    return(sml)
#}
#
#
#get.smile2 <- function(keggList) {
#    sml <- c()
#    for (i in 1:length(keggList)) {
#        keggId <- keggList[i] 
#        keggBget <- "http://www.genome.jp/dbget-bin/www_bget?" 
#        url <- paste(keggBget, keggId,  sep = "") 
#        txt <- getURL(url, write = basicTextGatherer()) 
#        v1 <- strsplit(txt, "\n")[[1]] 
#        p <- grep("sid=", v1, perl=TRUE) 
#        sid <- gsub("(^.+)sid=|\\\">\\d{1,10}.+$","",v1[p],perl=TRUE ) 
#        sidURL <- "http://pubchem.ncbi.nlm.nih.gov/summary/summary.cgi?sid="
#        url <- paste(sidURL, sid, "&viewopt=PubChem", sep = "") 
#        txt <- getURL(url, write = basicTextGatherer()) 
#        v2 <- grep("\\(CID \\d{1,10}\\)", strsplit(txt, "\n")[[1]], value=TRUE,perl=TRUE)
#        cid <- gsub("^.+ \\(CID |\\).+$","", v2,perl=TRUE) 
#        #p <- grep("chebi",v1, ignore.case=TRUE)[1]
#        #chebi <- gsub("^.+chebiId=CHEBI:|\\\">\\d{5}.+$", "", v1[p], perl=TRUE) 
#
#        cidURL <- "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?tool=rpubchem&db=pccompound&id=" 
#        #chebiURL <- "http://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI:"
#        url <- paste(cidURL, cid,  sep = "") 
#        txt <- getURL(url, write = basicTextGatherer()) 
#        v2 <- grep("CanonicalSmiles", strsplit(txt, "\n")[[1]], value=TRUE)
#        sml1 <- gsub("^.+String\\\">|</Item>$","",v2, perl=TRUE) 
#        if (length(sml1)) {
#            sml <- append(sml, sml1)
#         }
#        else {
#            sml <- append(sml, "miss")
#         }
#     }
#        return(sml)
#}
    if(length(descData==NULL)) {
#        if (get.smile==1) {
#             smlTest <- get.smile1(as.character(testSet$id))  
#             smlPred <- get.smile1(as.character(predSet$id))  
#        }
#        if (get.smile==2) {
#             registerDoMC() 
#             keggList <- as.character(testSet$id)
#             smlTest1 <- foreach(i=1:length(keggList), .combine=c) %dopar% {
#                            get.smile2(keggList[i])
#                         }
#             keggList <- predSet$id
#             smlPred1 <- foreach(i=1:length(keggList), .combine=c) %dopar% {
#                            get.smile2(keggList[i])
#                         }
#        }
# 
#        if (length(which(smlTest1=="miss"))) { smlTest <- smlTest1[-which(smlTest1=="miss")]; testSet <- testSet[-which(smlTest1=="miss"),] } 
#        if (length(which(smlPred1=="miss"))) { smlPred <- smlPred1[-which(smlPred1=="miss")]; predSet <- predSet[-which(smlModel1=="miss"),] }
#        if (sum(duplicated(smlTest))) { smlTest <- smlTest[-which(duplicated(smlTest))]; testSet <- testSet[-which(duplicated(smlPred)),] }
#        if (sum(duplicated(smlPred))) { smlPred <- smlPred[-which(duplicated(smlPred))]; predSet <- predSet[-which(duplicated(smlPred)),] }
#         molsTest <- parse.smiles(smlTest)
#         molsPred <- parse.smiles(smlPred)
#
#        if (is.null(descList)) {
#             none <- c("IPMolecularLearning", "AminoAcidCountDescriptor") 
#             selection <- get.desc.names()
#             c <- sapply(none, grep, selection) 
#             selection <- selection[-c] 
#             xT <- eval.desc(molsTest, selection, verbose=TRUE)
#             xP <- eval.desc(molsPred, selection, verbose=TRUE)
#        }
#        
#        else {
#            xT <- eval.desc(molsTest, descList, verbose=TRUE)
#            xP <- eval.desc(molsPred, descList, verbose=TRUE)
#        }
    }
    else {
        xT <-  descData$testM
        xP <-  descData$predM
    }
 modelFun <- function(x, testSet) {
       x2 <- x[,apply(x, 2, function(a) {all(!is.na(a))})]
       oc <- which(apply(x2, 2, sum)==0)
       if (length(oc)) x2 <- as.matrix(x2[, -oc])
	x3 <- x2
#       x3 <- as.matrix(apply(x2, 2, function(x) (x-mean(x))/sd(x)))
#       if (sum(is.na(x3))) x3 <- x3[,apply(x3, 2, function(a) {all(!is.na(a))})]
       coord1 <- c(); 
       for ( i in 1:(ncol(x3)-1) ) { 
          ld <- fixDependence(x3[,c(i,i+1)],x3[,-c(i,i+1)])
          if (length(!is.null(ld) & ld > i) ) { 
              coord1 <- append(coord1,ld) 
          }  
      }

      if (!is.null(coord1)) {  
          sel <- regsubsets(x=x3[,-unique(coord1)], y=testSet,nvmax=ncol(x3[,-unique(coord1)]),  method="exhaustive",really.big=T)
      } 
      else {  sel <- regsubsets(x=x3, y=testSet,nvmax=ncol(x3),  method="exhaustive",really.big=T)  
      }

     pred <- summary(sel) 
     pos <- which(pred$cp==min(pred$cp)) 
     n <- as.character(names(coef(sel, 1:length(pred$cp))[[pos]]))[-1]
     fm <- paste(n, collapse="+")  
     model <- lm(as.formula(paste("y~",fm, sep="")), data=data.frame(y=testSet,x3), y=TRUE )
     Y <- as.matrix(testSet)
     X <- as.matrix(x3[, n])
     theta.fit <- function(x,y){lsfit(x,y)} 
     theta.predict <- function(fit,x){cbind(1,x)%*%fit$coef}  
     results <- crossval(X,Y,theta.fit,theta.predict,ngroup=5)
     rawRLm <- cor(Y, model$fitted.values)**2 # raw R2 
     crossValRlm <- cor(Y,results$cv.fit)**2 # cross-validated R2 

        invisible(list(rawRLm=rawRLm, crossValRlm=crossValRlm, predError=model$residuals,coef=model$coefficients, fitted.values= model$fitted.values, n=n))
}
     # pred matrix does not need all descriptors
     myAns <- modelFun(xT, log10(testSet$data))
     #xP <- apply(xP, 2, function(x) (x-mean(x))/sd(x)) 
     pred <- as.matrix(cbind(1,xP[,myAns$n]),ncol=length(myAns$n)+1)%*%as.matrix(myAns$coef,ncol=1)

     invisible(list(ans=myAns, pred=pred, testSet=testSet, predSet=predSet, voidTime=voidTime))
}

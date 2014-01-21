#' sbml2table 
#'
#' Retrieves compound and its reactions from SBML models. 
#' @param file SBML file. 
#' @return A database in the format required for ProbMetab, with
#' mandatory fields id, formula and reactions. 
#' @export

sbml2table <-  
 function (file)  
 { 
#     require(XML) 
     doc <- xmlParse(file, ignoreBlanks = TRUE) 
     name <- xpathSApply(doc, "//name:species", xmlGetAttr, "name",  
         namespaces = "name") 
     formula <- xpathSApply(doc, "//formula:species", xmlGetAttr,  
         "formula", namespaces = "formula") 
     id.sbml <- xpathSApply(doc, "//id:species", xmlGetAttr, "id",  
         namespaces = "id") 
     reacs <- xpathSApply(doc, "//id:reaction", xmlGetAttr, "id",  
         namespaces = "id") 
#     nodesReactants <- getNodeSet(doc, "//speciesReference:listOfReactants",  
#         namespaces = "speciesReference") 
#     lreactants <- sapply(nodesReactants, function(x) sapply(xmlChildren(x),  
#         xmlGetAttr, "species")) 
#     nodesProducts <- getNodeSet(doc, "//speciesReference:listOfProducts",  
#         namespaces = "speciesReference") 
#     lproducts <- sapply(nodesProducts, function(x) sapply(xmlChildren(x),  
#         xmlGetAttr, "species")) 
#     allreacs <- c("", "") 
#     for (i in 1:length(reacs)) { 
#         allreacs <- rbind(allreacs, cbind(reacs[i], lreactants[[i]]),  
#             cbind(reacs[i], lproducts[[i]])) 
#     } 
#
Reacts <- getNodeSet(doc, "//speciesReference:reaction",  namespaces = "speciesReference") 
     allreacs <- c("", "") 
     for (i in 1:length(reacs)) { 

lreac <-""
lprod <-"" 
child <- xmlChildren(Reacts[[i]])

if(!is.null(child$listOfReactants)){
lreac <- unlist(lapply(lapply(as.list(xmlChildren(child$listOfReactants)), xmlAttrs), function(x) x["species"]) ) 
}
if(!is.null(child$listOfProducts)){
lprod <- unlist(lapply(lapply(as.list(xmlChildren(child$listOfProducts)), xmlAttrs), function(x) x["species"]) ) 
}

         allreacs <- rbind(allreacs, cbind(reacs[i], lreac),  
             cbind(reacs[i], lprod)) 
}
allreacs <- allreacs[-which(allreacs[,1]=="" | allreacs[,2]==""),] 
     allreacs <- unique(allreacs) 
     unqreac <- paste("r", sprintf("%04d", 1:length(unique(allreacs[,  
         1]))), sep = "") 
     allreacs2 <- rep(unqreac, table(allreacs[, 1])[unique(allreacs[,  
         1])]) 
     reactions.sbml <- sapply(id.sbml, function(x) paste(allreacs[which(allreacs[,  
         2] == x), 1], collapse = ";")) 
     reactions <- sapply(id.sbml, function(x) paste(allreacs2[which(allreacs[,  
         2] == x)], collapse = ";")) 
     cpdData <- cbind(id.sbml, name, formula, reactions, reactions.sbml) 
     cpdData[, 1] <- sub("(.+)_.+$", "\\1", cpdData[, 1]) 
     cpdData <- unique(cpdData) 
     id <- paste("c", sprintf("%04d", 1:length(cpdData[, 1])),  
         sep = "") 
     cpdData <- cbind(id, cpdData) 
     rownames(cpdData) <- 1:nrow(cpdData) 
     return(cpdData) 
     # using rsbml 
 #    require(rsbml) 
 #    sbml <- rsbml_read(file) 
 #    v <- rep("", 3) 
 #    for (i in 1:length(sbml@model@species)) { 
 #        id <- sbml@model@species[[i]]@id 
 #        nf1 <- strsplit(sbml@model@species[[i]]@name, "_")[[1]] 
 #        nf <- nf1[nf1 != ""] 
 #        v <- rbind(v, c(id, nf)) 
 #    } 
 #    v <- v[-1, ] 
 #    v[, 1] <- sub("_\\w$", "", v[, 1]) 
 #    v <- unique(v) 
 #    ln <- sapply(v[, 1], grep, v[, 1]) 
 #    if (sum(unlist(lapply(ln, length)) > 1))  
 #        cat("\nWarning: non-unique ids\n") 
 #    r <- rep("", 2) 
 #    for (i in 1:length(sbml@model@reactions)) { 
 #        cpds <- unique(c(names(sbml@model@reactions[[i]]@reactants),  
 #            names(sbml@model@reactions[[i]]@products))) 
 #        reac <- cbind(cpds, sbml@model@reactions[[i]]@id) 
 #        r <- rbind(r, reac) 
 }

#' formula2mass 
#'
#' Calculates the exact mass of a given formula. 
#' @param formula molecular formula, not allowing characters aside
#' atom's alphabet.
#' @return a float representing the sum of monoisotopic masses. 
#'
#' @export

formula2mass <-
function(formula) {
    # very simple function, avoid using that
    f.extract <- function(formula) 
    { 
         # pattern to match the initial chemical 
         # assumes chemical starts with an upper case and optional lower case followed 
         # by zero or more digits. 
         first <- "^([[:upper:]][[:lower:]]?)([0-9]*).*" 
         # inverse of above to remove the initial chemical 
         last <- "^[[:upper:]][[:lower:]]?[0-9]*(.*)" 
         result <- list() 
         extract <- formula 
         # repeat as long as there is data 
         while ((start <- nchar(extract)) > 0){ 
             chem <- sub(first, '\\1 \\2', extract) 
             extract <- sub(last, '\\1', extract) 
             # if the number of characters is the same, then there was an error 
             if (nchar(extract) == start){ 
                 warning("Invalid formula:", formula) 
                 return(NULL) 
             } 
             # append to the list 
             result[[length(result) + 1L]] <- strsplit(chem, ' ')[[1]] 
         } 
         result 
     } 
    #data <- system.file("data/NIST_relative_atomic_mass.rda", package="probmetab")
    #load(data) 
    atom <- NIST_relative_atomic_mass
    byMass <- f.extract(formula) 
    if(is.null(byMass)) return(NULL)
    key <- 0
    lapply(byMass, function(x)  if(!as.logical(length(which(atom[,1]==x[1])))) key <<- 1)
    if (key==0) {
        listOfMasses <- lapply(byMass, function(x) 
                                       if(!is.na(x[2])) atom[which(atom[,1]==x[1]),2]*as.numeric(x[2]) 
                                       else atom[which(atom[,1]==x[1]),2]
                               )

        return(sum(unlist(listOfMasses)))
    }
    else {
          warning("Invalid formula:", formula) 
          return(NULL) 
    }
}

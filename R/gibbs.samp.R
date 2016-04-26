#' gibbs.samp
#'
#' Call a C++ function that runs a Gibbs Sampler algorithm  
#' to sample from the distribution of metabolite to compound attribution
#' with the previous assumption that the connected combination of
#' attributions makes more sense. 
#'
#' @param x a vector of masses (unique from mass/retention time pairs).
#' @param y a vector of candidate compounds for each mass. 
#' @param N number of iterations to sample.
#' @param w matrix of compound connections.
#' @param p matrix of likelihood probabilities.
#' @return A list of matrices with attribution indexes and probabilities.
#'
#' @export 
#' 
 gibbs.samp <- function (x, y, N, w, p) 
.Call("file193b1b67af14", x, y, N, w, p, PACKAGE = "ProbMetab")


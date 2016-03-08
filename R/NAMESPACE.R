#' @importFrom assertthat assert_that
#' @importFrom magrittr %>% %<>%
#' @importFrom stats weights
#' @import aoos
#' @import methods
#' @import Matrix
#' @import modules
#' @import Rcpp
#' @useDynLib saeRobust
NULL

globalVariables(c(".", ".self", "psi", "convCrit", "maxIter", "maxIterParam", "k", "K"))

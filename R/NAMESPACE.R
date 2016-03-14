#' @importFrom assertthat assert_that
#' @importFrom magrittr %>% %<>%
#' @importFrom stats weights getCall formula
#' @import Matrix
#' @import aoos
#' @import methods
#' @import modules
#' @import Rcpp
#' @useDynLib saeRobust
NULL

globalVariables(c(".", ".self", "psi", "convCrit", "maxIter", "maxIterParam", "k", "K"))

# needed for update because it is an S4 method in Matrix:
setOldClass("rfh")
setOldClass(c("fitrsfh", "fitrfh"))
setOldClass(c("fitrtfh", "fitrfh"))
setOldClass(c("fitrstfh", "fitrfh"))

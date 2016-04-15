#' @importFrom assertthat assert_that
#' @importFrom magrittr %>% %<>%
#' @importFrom stats weights getCall formula
#' @import Matrix
#' @import aoos
#' @import methods
#' @import modules
#' @import Rcpp
#' @import ggplot2
#' @useDynLib saeRobust
NULL

globalVariables(c(".", ".self", "psi", "convCrit", "maxIter", "maxIterParam", "k", "K"))

# needed for S4-dispatch:
setOldClass("rfh")
setOldClass("rfhVariance")
setOldClass(c("fitrsfh", "fitrfh"))
setOldClass(c("fitrtfh", "fitrfh"))
setOldClass(c("fitrstfh", "fitrfh"))

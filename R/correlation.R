#' Correlation Structure
#'
#' Various correlation structures. They can be used inside the \link{rfh}
#' function to supply an alterantive variance structure to be fitted. For
#' examples see the documentation of \link{rfh}.
#'
#' @details \code{corSAR1} can be used to model a simultanous autoregressive
#'   process of order one: spatial correlation.
#'
#' @slot W the row-standardised proximity matrix
#'
#' @rdname correlation
#' @export
setClass("corSAR1", slots = c(W = "matrixORMatrix"))

#' @name corSAR1
#' @usage corSAR1(W)
#'
#' @param W the row-standardised proximity matrix
#'
#' @rdname correlation
#' @export corSAR1
corSAR1 <- function(W) {
  new("corSAR1", W = W)
}

#' @details \code{corAR1} can be used to model a autoregressive
#'   process of order one: temporal correlation.
#'
#' @slot nTime (numeric) number of time periods
#'
#' @rdname correlation
#' @export
setClass("corAR1", slots = c(nTime = "numericORinteger"))

#' @name corAR1
#' @usage corAR1(nTime)
#'
#' @param nTime (numeric) number of time periods
#'
#' @rdname correlation
#' @export corAR1
corAR1 <- function(nTime) {
  new("corAR1", nTime = nTime)
}

#' @rdname correlation
#' @export
setClass("corSAR1AR1", contains = c("corAR1", "corSAR1"))

#' @details \code{corSAR1AR1} can be used to model to model spatial and temporal
#'   correlation
#'
#' @name corSAR1AR1
#' @usage corSAR1AR1(nTime, W)
#'
#' @rdname correlation
#' @export corSAR1AR1
corSAR1AR1 <- function(nTime, W) {
  new("corSAR1AR1", nTime = nTime, W = W)
}

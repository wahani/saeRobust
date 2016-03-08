#' Correlation Structure
#'
#' @slot W the proximity matrix
#'
#' @rdname correlation
#' @export
corSAR1(W ~ matrix | Matrix) %type% {
  .Object
}

#' @name corSAR1
#' @usage corSAR1(W, ...)
#'
#' @param W the proximity matrix
#' @param ... arguments passed to \code{new}
#'
#' @rdname correlation
#' @export corSAR1
corSAR1

#' @slot nTime (numeric) number of time periods
#'
#' @rdname correlation
#' @export
corAR1(nTime ~ numeric | integer) %type% {
  .Object
}

#' @name corAR1
#' @usage corAR1(nTime, ...)
#'
#' @param nTime (numeric) number of time periods
#'
#' @rdname correlation
#' @export corAR1
corAR1

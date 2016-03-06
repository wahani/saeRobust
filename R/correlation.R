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

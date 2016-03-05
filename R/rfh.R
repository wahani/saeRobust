#' Robust Fay Herriot Model
#'
#' @param formula (formula)
#' @param data (data.frame)
#' @param samplingVar (character)
#' @param ... arguments passed to methods
#' @param x0 (numeric) starting values for variance parameters
#' @param k (numeric) tuning constant
#' @param tol (numeric) numerical toloerance to be used during optimisation
#' @param y (numeric) response vector
#' @param x ([m|M]atrix) the design matrix
#' @param maxIter (integer) the maximum number of iterations
#' @param maxIterRe (integer) the maximum number of iterations for fitting the
#'   random effects
#'
#' @rdname rfh
#'
#' @export
rfh(formula, data, samplingVar, ...) %g% standardGeneric("rfh")

#' @rdname rfh
#' @export
rfh(formula ~ formula, data ~ data.frame, samplingVar ~ character, ...) %m% {
  call <- match.call()
  xy <- makeXY(formula, data)
  samplingVar <- data[[samplingVar]]

  retList(
    "rfh",
    public = c("call"),
    super = fitrfh(xy$y, xy$x, samplingVar, ...)
  )

}

#' @param object (rfh) an object of class rfh
#' @param type (character) one in \code{c("linear", "REBLUP")}
#'
#' @rdname rfh
#' @export
predict.fitrfh <- function(object, type = "REBLUP", ...) {

  re <- if ("linear" == type) {
    0
  } else if ("REBLUP" == type) {
    as.numeric(variance(object)$Z() %*% object$re)
  }

  Xb <- object$x %*% object$coefficients
  out <- data.frame(prediction = as.numeric(Xb + re))
  names(out) <- type
  out$re <- re

  out

}

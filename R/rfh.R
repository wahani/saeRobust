#' Robust Fay Herriot Model
#'
#' @param formula (formula)
#' @param data (data.frame)
#' @param samplingVar (character)
#' @param correlation an optional correlation structure, e.g. \link{corSAR1}.
#' @param ... arguments passed \link{fitGenericModel}
#'
#' @rdname rfh
#'
#' @export
rfh(formula, data, samplingVar, correlation = NULL, ...) %g% standardGeneric("rfh")

#' @name rfh
#' @usage \S4method{rfh}{formula,data.frame,character,ANY}(formula, data, samplingVar, correlation, ...)
#' @aliases rfh,formula,data.frame,character,ANY-method
#' @rdname rfh
NULL

rfh(formula ~ formula, data ~ data.frame, samplingVar ~ character, correlation ~ ANY, ...) %m% {
  call <- match.call()
  xy <- makeXY(formula, data)
  samplingVar <- data[[samplingVar]]

  stripSelf(retList(
    "rfh",
    public = c("call"),
    super = rfh(xy$y, xy$x, samplingVar, correlation, ...)
  ))

}

#' @rdname rfh
#' @export
rfh(formula ~ numeric, data ~ matrix | Matrix, samplingVar ~ numeric, correlation ~ 'NULL', ...) %m% {
  fitrfh(formula, data, samplingVar, ...)
}

#' @rdname rfh
#' @export
rfh(formula ~ numeric, data ~ matrix | Matrix, samplingVar ~ numeric, correlation ~ corSAR1, ...) %m% {
  fitrsfh(formula, data, samplingVar, correlation@W, ...)
}

#' @rdname rfh
#' @export
rfh(formula ~ numeric, data ~ matrix | Matrix, samplingVar ~ numeric, correlation ~ corAR1, ...) %m% {
  fitrtfh(formula, data, samplingVar, correlation@nTime, ...)
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

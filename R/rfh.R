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
    public = c("call", "formula"),
    super = rfh(xy$y, xy$x, samplingVar, correlation, ...)
  ))

}

#' @rdname fit
#' @export
rfh(formula ~ numeric, data ~ matrix | Matrix, samplingVar ~ numeric, correlation ~ 'NULL', ...) %m% {
  fitrfh(formula, data, samplingVar, ...)
}

#' @rdname fit
#' @export
rfh(formula ~ numeric, data ~ matrix | Matrix, samplingVar ~ numeric, correlation ~ corSAR1, ...) %m% {
  fitrsfh(formula, data, samplingVar, correlation@W, ...)
}

#' @rdname fit
#' @export
rfh(formula ~ numeric, data ~ matrix | Matrix, samplingVar ~ numeric, correlation ~ corAR1, ...) %m% {
  fitrtfh(formula, data, samplingVar, correlation@nTime, ...)
}

#' @rdname fit
#' @export
rfh(formula ~ numeric, data ~ matrix | Matrix, samplingVar ~ numeric, correlation ~ corSAR1AR1, ...) %m% {
  fitrstfh(formula, data, samplingVar, correlation@W, correlation@nTime, ...)
}

#' @param object (rfh) an object of class rfh
#' @param type (character) one or more in \code{c("linear", "reblup")}
#'
#' @rdname rfh
#' @export
predict.fitrfh <- function(object, type = "reblup", ...) {

  re <- as.numeric(variance(object)$Z() %*% object$re)
  Xb <- fitted.values(object)

  out <- data.frame(re = re)
  if (is.element("linear", type)) out$linear <- Xb
  if (is.element("reblup", type)) out$reblup <- Xb + re

  out

}

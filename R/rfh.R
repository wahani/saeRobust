#' Robust Fay Herriot Model
#'
#' @param formula (formula)
#' @param data (data.frame)
#' @param samplingVar (character)
#' @param ... arguments passed \link{fitGenericModel}
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

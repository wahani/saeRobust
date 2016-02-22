#' Construct variance
#'
#' A generic function to construct the different variance components of an
#' object.
#'
#' @param .object,object an object
#' @param re (numeric | NULL) if NULL then the random effects have to be
#'   computet.
#' @param ... arguments passed to method
#'
#' @export
#' @rdname variance
variance <- function(.object, ...) UseMethod("variance")

#' @export
#' @rdname variance
variance.rfh <- function(.object, ...) {

    expose(matVFH(.object$variance, .object$samplingVar))

    retList("rfhVariance")
}

#' @export
#' @rdname variance
weights.rfh <- function(object, re = NULL, ...) {

    V <- variance(object)
    re <- if (is.null(re)) predict(object)[["re"]] else re

    W <-  matW(
        y = object$xy$y,
        X = object$xy$x,
        beta = object$beta,
        u = re,
        matV = V,
        psi = object$psi
    )

    A <- matA(
        y = object$xy$y,
        X = object$xy$x,
        beta = object$beta,
        V = V$V(),
        Vinv = V$VInv(),
        psi = object$psi
    )

    B <- matB(
        y = object$xy$y,
        X = object$xy$x,
        beta = object$beta,
        u = re,
        V,
        psi = object$psi
    )

    retList("rfhWeights", c("W", "A", "B")) %>% stripSelf

}


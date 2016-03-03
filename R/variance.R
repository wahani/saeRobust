#' Construct variance
#'
#' A generic function to construct the different variance components of an
#' object.
#'
#' @param .object,object an object
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
weights.rfh <- function(object, ...) {

    V <- variance(object)

    W <-  matW(
        y = object$y,
        X = object$x,
        beta = object$coefficients,
        u = object$re,
        matV = V,
        psi = object$psi
    )

    A <- matA(
        y = object$y,
        X = object$x,
        beta = object$coefficients,
        matV = V,
        psi = object$psi
    )

    B <- matB(
        y = object$y,
        X = object$x,
        beta = object$coefficients,
        u = object$re,
        V,
        psi = object$psi
    )

    stripSelf(retList("rfhWeights", c("W", "A", "B")))

}


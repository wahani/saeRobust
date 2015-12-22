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

    V <- function() {
        Diagonal(x = .object$samplingVar + .object$variance)
    }

    Vinv <- function() {
        solve(.self$V())
    }

    Vu <- function() {
        Diagonal(length(.object$samplingVar), .object$variance)
    }

    Ve <- function() {
        Diagonal(x = .object$samplingVar)
    }

    VuInv <- function() {
        solve(.self$Vu())
    }

    VeInv <- function() {
        solve(.self$Ve())
    }

    VuSqrtInv <- function() {
        sqrt(.self$VuInv())
    }

    VeSqrtInv <- function() {
        sqrt(.self$VeInv())
    }

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
        V = V$V(),
        Vinv = V$Vinv(),
        VuSqrtInv = V$VuSqrtInv(),
        VeSqrtInv = V$VeSqrtInv(),
        psi = object$psi
    )

    A <- matA(
        y = object$xy$y,
        X = object$xy$x,
        beta = object$beta,
        V = V$V(),
        Vinv = V$Vinv(),
        psi = object$psi
    )

    B <- matB(
        y = object$xy$y,
        X = object$xy$x,
        beta = object$beta,
        u = re,
        VuSqrtInv = V$VuSqrtInv(),
        VeSqrtInv = V$VeSqrtInv(),
        psi = object$psi
    )

    retList("rfhWeights", c("W", "A", "B")) %>% stripSelf

}


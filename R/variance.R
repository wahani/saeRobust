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
variance.fitrfh <- function(.object, ...) {

    expose(matVFH(.object$variance, .object$samplingVar))

    retList("rfhVariance")
}

#' @export
#' @rdname variance
variance.fitrsfh <- function(.object, ...) {

  expose(matVSFH(
    .object$variance[1],
    .object$variance[2],
    .object$W,
    .object$samplingVar))

  retList("rfhVariance")
}

#' @export
#' @rdname variance
variance.fitrtfh <- function(.object, ...) {

  expose(matVTFH(
    .object$variance[1],
    .object$variance[c(2, 3)],
    .object$nTime,
    .object$samplingVar))

  retList("rfhVariance")
}

#' @export
#' @rdname variance
variance.fitrstfh <- function(.object, ...) {

  expose(matVSTFH(
    .object$variance[1:2],
    .object$variance[3:4],
    .object$W,
    .object$nTime,
    .object$samplingVar))

  retList("rfhVariance")

}

#' @export
#' @rdname variance
weights.fitrfh <- function(object, ...) {

    V <- variance(object)

    W <-  matW(
        y = object$y,
        X = object$x,
        beta = object$coefficients,
        u = object$re,
        matV = V,
        psi = . %>% psiOne(object$k)
    )

    A <- matA(
        y = object$y,
        X = object$x,
        beta = object$coefficients,
        matV = V,
        psi = . %>% psiOne(object$k)
    )

    B <- matB(
        y = object$y,
        X = object$x,
        beta = object$coefficients,
        u = object$re,
        V,
        psi = . %>% psiOne(object$k)
    )

    stripSelf(retList("rfhWeights", c("W", "A", "B")))

}


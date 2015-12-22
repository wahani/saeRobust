#' Compute the Mean Squared Error of an Estimator
#'
#' A generic function to compute the mean squared error of the predicted values
#' under the estimated model. This estimation can also be triggered using the
#' predict function.
#'
#' @param object (see methods) an object containing the estimation result, e.g.
#'   \link{rfh}
#' @param type (character) the type of the MSE. Available are 'pseudo' and
#'   'boot'.
#' @param re (NULL | numeric) optionally provide the estimated random effects
#' @param V (NULL | list) optionally provide the variance using \link{variance}
#' @param ... arguments passed to methods
#'
#' @export
#' @rdname mse
mse <- function(object, ...) UseMethod("mse")

#' @export
#' @rdname mse
mse.rfh <- function(object, type = "pseudo", re = NULL, V = NULL, ...) {

    pseudo <- function() {
        W <- weights(object, re)$W
        as.numeric(W^2 %*% diag(V$V()) + (W %*% Xb - Xb)^2)
    }

    boot <- function() {
        # TODO: implement bootstrap MSE
        NULL
    }

    re <- if (is.null(re)) predict(object)[["re"]] else re
    V <- if (is.null(V)) variance(object) else V
    Xb <- object$xy$x %*% object$beta

    out <- data.frame(REBLUP = as.numeric(Xb + re))
    if ("pseudo" %in% type) out$pseudo <- pseudo()
    if ("boot" %in% type) out$boot <- boot()
    out

}

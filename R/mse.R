#' Compute the Mean Squared Error of an Estimator
#'
#' A generic function to compute the mean squared error of the predicted values
#' under the estimated model. This estimation can also be triggered using the
#' predict function.
#'
#' @param object (see methods) an object containing the estimation result, e.g.
#'   \link{rfh}
#' @param type (character) the type of the MSE. Available are 'pseudo' and
#'   'boot'
#' @param re (NULL | numeric) optionally provide the estimated random effects
#' @param V (NULL | list) optionally provide the variance using \link{variance}
#' @param B (numeric) number of bootstrap repetitions
#' @param ... arguments passed to methods
#'
#' @export
#' @rdname mse
mse <- function(object, ...) UseMethod("mse")

#' @export
#' @rdname mse
mse.rfh <- function(object, type = "pseudo", re = NULL, V = NULL, B = 100, ...) {

    pseudo <- function() {
        W <- weights(object, re)$W
        as.numeric(W^2 %*% diag(V$V()) + (W %*% Xb - Xb)^2)
    }

    boot <- function(B) {
        bootSamples <- replicate(B, {
            re <- MASS::mvrnorm(1, mu = rep(0, length(Xb)), V$Vu())
            e <- rnorm(length(Xb), mean = 0, sd = sqrt(object$samplingVar))
            trueVal <- Xb + re
            y <- trueVal + e
            fit <- fitrfh(
                y = y,
                X = object$xy$x,
                samplingVar = object$samplingVar,
                theta0 = c(object$beta, object$variance),
                psi = object$psi
            )
            (predict(fit)$REBLUP - trueVal)^2
        })
        rowMeans(bootSamples)
    }

    re <- if (is.null(re)) predict(object)$re else re
    V <- if (is.null(V)) variance(object) else V
    Xb <- as.numeric(object$xy$x %*% object$beta)

    out <- data.frame(REBLUP = Xb + re)
    if ("pseudo" %in% type) out$pseudo <- pseudo()
    if ("boot" %in% type) out$boot <- boot(B)
    out

}

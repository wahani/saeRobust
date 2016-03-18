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
#' @param B (numeric) number of bootstrap repetitions
#' @param ... arguments passed to methods
#'
#' @export
#' @rdname mse
mse <- function(object, ...) UseMethod("mse")

#' @export
#' @rdname mse
mse.fitrfh <- function(object, type = "pseudo", B = 100, ...) {

    pseudo <- function(Xb, matV) {
        W <- weights(object)$W
        as.numeric(W^2 %*% diag(matV$V()) + (W %*% Xb - Xb)^2)
    }

    boot <- function(B, Xb, matV) {
        bootSamples <- replicate(B, {
            re <- MASS::mvrnorm(1, mu = rep(0, length(Xb)), matV$Vu())
            e <- MASS::mvrnorm(1, mu = rep(0, length(Xb)), matV$Ve())
            trueVal <- Xb + matV$Z() %*% re
            y <- trueVal + e
            fit <- fitrfh(
                y = y,
                x = object$x,
                samplingVar = object$samplingVar,
                x0Coef = object$coefficients,
                x0Var = object$variance,
                k = object$k,
                tol = object$tol
            )
            (fit$reblup - as.numeric(trueVal))^2
        })
        rowMeans(bootSamples)
    }

    matV <- variance(object)
    Xb <- fitted.values(object)

    out <- predict(object)["reblup"]
    if ("pseudo" %in% type) out$pseudo <- pseudo(Xb, matV)
    if ("boot" %in% type) out$boot <- boot(B, Xb, matV)
    out

}

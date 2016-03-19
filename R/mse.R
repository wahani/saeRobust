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
    G <- diag(matV$Z() %*% tcrossprod(matV$Vu(), matV$Z()))
    A <- (W - Diagonal(length(object$y)))^2
    var <- A %*% G + W^2 %*% object$samplingVar
    bias <- W %*% Xb - Xb
    as.numeric(var + bias^2)
  }

  boot <- function(B, matV) {
    resList <- bootstrap(object, matV, B = B, filter = c("reblup", "trueY"))
    colMeans(do.call(rbind, lapply(resList, function(el) {
      (el$reblup - el$trueY)^2
    })))
  }

  matV <- variance(object)
  Xb <- fitted.values(object)

  out <- predict(object)["reblup"]
  if ("pseudo" %in% type) out$pseudo <- pseudo(Xb, matV)
  if ("boot" %in% type) out$boot <- boot(B, matV)
  out

}

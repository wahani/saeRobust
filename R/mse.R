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
#' @param predType (character) the type of prediction: \code{c("reblup", "reblupbc")}
#' @param ... arguments passed to methods
#'
#' @export
#' @rdname mse
mse <- function(object, ...) UseMethod("mse")

#' @export
#' @rdname mse
mse.fitrfh <- function(object, type = "pseudo", predType = "reblupbc", B = 100, ...) {

  pseudo <- function(Xb, matV, Wtype = "W") {
    W <- weights(object)[[Wtype]]
    G <- diag(matV$Z() %*% tcrossprod(matV$Vu(), matV$Z()))
    A <- (W - Diagonal(length(object$y)))^2
    var <- A %*% G + W^2 %*% object$samplingVar
    bias <- if (Wtype == "W") W %*% Xb - Xb else 0
    as.numeric(var + bias^2)
  }

  boot <- function(B, matV) {

    resList <- bootstrap(
      object, matV, B = B,
      postProcessing = function(fit) c(predict(fit, type = c("reblup", "reblupbc")), fit["trueY"])
    )

    reblup <- colMeans(do.call(rbind, lapply(resList, function(el) {
      (el$reblup - el$trueY)^2
    })))

    reblupbc <- colMeans(do.call(rbind, lapply(resList, function(el) {
      (el$reblupbc - el$trueY)^2
    })))

    retList(public = c("reblup", "reblupbc"))

  }

  matV <- variance(object)
  Xb <- fitted.values(object)

  out <- predict(object, type = predType)[predType]
  if ("pseudo" %in% type && "reblup" %in% predType) out$pseudo <- pseudo(Xb, matV)
  if ("pseudo" %in% type && "reblupbc" %in% predType) out$pseudobc <- pseudo(Xb, matV, "Wbc")
  if ("boot" %in% type) {
    boots <- boot(B, matV)
    if ("reblup" %in% predType) out$boot <- boots$reblup
    if ("reblupbc" %in% predType) out$bootbc <- boots$reblupbc
  }
  out

}

#' Fit model on Bootstrap sample
#'
#' These methods help to repeatedly fit a \link{rfh} model on bootstrap samples.
#'
#' @param modelA a fitted object
#' @param modelB a fitted object used to draw samples. In most cases this is
#'   \code{modelA}. Alternatively it may be useful to use a non-robust model.
#' @param B the number of repetitions
#' @param filter a vector indicating which elements in the fittedd object to
#'   keep in each repetition.
#' @param ... arguments passed down to methods
#'
#' @export
#' @rdname bootstrap
bootstrap(modelA, modelB = modelA, B = NULL, ...) %g% standardGeneric("bootstrap")

#' @export
#' @rdname bootstrap
bootstrap(modelA, modelB ~ missing, B ~ missing, ...) %m%
  bootstrap(modelA = modelA, modelB = modelA, B = NULL, ...)

#' @export
#' @rdname bootstrap
bootstrap(modelA, modelB ~ missing, B ~ numeric|integer, ...) %m% {
  bootstrap(modelA = modelA, modelB = modelA, B = B, ...)
}

#' @export
#' @rdname bootstrap
bootstrap(modelA, modelB, B ~ integer|numeric, filter = NULL, ...) %m% {
  if (is.null(filter)) {
    replicate(B, bootstrap(modelA, modelB, NULL, ...), FALSE)
  } else {
    replicate(B, bootstrap(modelA, modelB, NULL, ...)[filter], FALSE)
  }
}

#' @export
#' @rdname bootstrap
bootstrap(modelA ~ rfh, modelB ~ rfh, B ~ NULL, ...) %m% {
  # This we need to get directly in the update method for fitrfh class
  # Otherwise weired things are happening in the call for the S4-generic
  class(modelA) <- class(modelA)[-1] # this hopefully only removes 'rfh'

  # Bootstrap sample:
  Xb <- fitted.values(modelB)
  matV <- variance(modelB)
  re <- MASS::mvrnorm(1, mu = rep(0, length(Xb)), matV$Vu())
  e <- MASS::mvrnorm(1, mu = rep(0, length(Xb)), matV$Ve())
  trueY <- as.numeric(Xb + matV$Z() %*% re)
  y <- trueY + e

  # refit:
  out <- update(modelA, y = y)
  out$trueY <- trueY
  out
}


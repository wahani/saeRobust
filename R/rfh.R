#' Robust Fay Herriot Model
#'
#' @param formula (formula)
#' @param data (data.frame)
#' @param samplingVar (character)
#' @param ... arguments passed to methods
#' @param x0 (numeric) starting values for variance parameters
#' @param k (numeric) tuning constant
#' @param tol (numeric) numerical toloerance to be used during optimisation
#' @param y (numeric) response vector
#' @param x ([m|M]atrix) the design matrix
#' @param maxIter (integer) the maximum number of iterations
#'
#' @rdname rfh
#'
#' @export
rfh(formula, data, samplingVar, ...) %g% standardGeneric("rfh")

#' @rdname rfh
#' @export
rfh(formula ~ formula, data ~ data.frame, samplingVar ~ character, ...) %m% {
  call <- match.call()
  xy <- makeXY(formula, data)
  samplingVar <- data[[samplingVar]]

  retList(
    public = c("call"),
    super = fitrfh(xy$y, xy$x, samplingVar, ...)
  )

}

#' @rdname rfh
#' @export
fitrfh <- function(y, x, samplingVar, x0 = 1, k = 1.345, tol = 1e-6, maxIter = 100) {
  # Non interactive fitting function for robust FH
  # y: (numeric) response
  # x: ((M|m)atrix) design matrix
  # samplingVar: (numeric) the sampling variances
  # x0: (numeric) starting values
  # k: (numeric) tuning constant

  # robust settings:
  psi <- . %>% psiOne(k)
  K <- getK(k)

  # Algorithm settings
  convCrit <- convCritRelative(tol)

  oneIter <- function(param) {
    param <- as.numeric(param)
    beta <- param[-length(param)]
    sigma2 <- param[length(param)]

    fpBeta <- fixedPointRobustBeta(
      y, x, matVFH(sigma2, samplingVar), psi = psi
    ) %>% addHistory

    beta <- fixedPoint(fpBeta, beta, addMaxIter(convCrit, maxIter))

    fpSigma2 <- fixedPointRobustDelta(
      y, x, beta, . %>% matVFH(samplingVar), psi, K
    ) %>% addConstraintMin(0) %>% addHistory

    sigma2 <- fixedPoint(fpSigma2, sigma2, addMaxIter(convCrit, maxIter))

    list(beta, sigma2)
  }

  # Fitting Model Parameter:
  x0 <- c(lm.fit(as.matrix(x), as.numeric(y))$coefficients, x0)
  out <- fixedPoint(addStorage(oneIter), x0, addMaxIter(convCrit, maxIter))
  iterations <- storage$reformat(attr(out, "storage"))
  coefficients <- out[-length(out)]
  names(coefficients) <- colnames(x)
  variance <- out[length(out)]
  names(variance) <- "var"

  # Fitting Random Effects
  re <- fitRe(y, x, coefficients, matVFH(variance, samplingVar), psi, convCrit)
  iterations <- c(iterations, list(attr(re, "history")))
  re <- as.numeric(re)
  names(iterations) <- c("coefficients", "variance", "re")

  stripSelf(retList("rfh", public = c(
    "coefficients", "variance", "psi", "samplingVar", "y", "x", "iterations",
    "k", "tol", "K", "re"
  )))

}

#' @param object (rfh) an object of class rfh
#' @param type (character) one in \code{c("linear", "REBLUP")}
#'
#' @rdname rfh
#' @export
predict.rfh <- function(object, type = "REBLUP", ...) {

  re <- if ("linear" == type) {
    0
  } else if ("REBLUP" == type) {
    as.numeric(variance(object)$Z() %*% object$re)
  }

  Xb <- object$x %*% object$coefficients
  out <- data.frame(prediction = as.numeric(Xb + re))
  names(out) <- type
  out$re <- re

  out

}

fitRe <- function(y, x, beta, matV, psi, convCrit) {
  # y: (numeric) response
  # x: (Matrix) Design-Matrix
  # beta: (numeric) estimated fixed effects
  # matV: (list) list of functions see e.g. matVFH
  # psi: (function) influence function
  # convCrit: (function) convergence criterion

  # Non-robust random effects as starting values:
  startingValues <- as.numeric(
    tcrossprod(matV$Vu(), matV$Z()) %*% matV$VInv() %*% (y - x %*% beta)
  )

  fpFun <- fixedPointRobustRandomEffect(
    y, x, beta, matV, psi
  )

  fixedPoint(addHistory(fpFun), startingValues, convCrit)

}

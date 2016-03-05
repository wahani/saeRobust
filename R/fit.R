#' @rdname rfh
#' @export
fitrfh <- function(
  y, x, samplingVar,
  k = 1.345, K = getK(k), psi = . %>% psiOne(k),
  x0Coef = NULL, x0Var = 1, x0Re = NULL,
  tol = 1e-6, maxIter = 100, maxIterRe = 100, convCrit = convCritRelative(tol)) {
  # Non interactive fitting function for robust FH
  # y: (numeric) response
  # x: ((M|m)atrix) design matrix
  # samplingVar: (numeric) the sampling variances
  # x0: (numeric) starting values
  # k: (numeric) tuning constant

  oneIter <- function(param) {
    param <- as.numeric(param)
    beta <- param[-length(param)]
    sigma2 <- param[length(param)]

    fpBeta <- addHistory(fixedPointRobustBeta(
      y, x, matVFH(sigma2, samplingVar), psi = psi
    ))

    beta <- fixedPoint(fpBeta, beta, addMaxIter(convCrit, maxIter))

    fpSigma2 <- fixedPointRobustDelta(
      y, x, beta, . %>% matVFH(samplingVar), psi, K
    ) %>% addConstraintMin(0) %>% addHistory

    sigma2 <- fixedPoint(fpSigma2, sigma2, addMaxIter(convCrit, maxIter))

    list(beta, sigma2)
  }

  # Fitting Model Parameter:
  if (is.null(x0Coef)) {
    x0Coef <- as.numeric(fitCoefStartingValue(y, x, matVFH(x0Var, samplingVar)))
  }
  out <- fixedPoint(addStorage(oneIter), c(x0Coef, x0Var), addMaxIter(convCrit, maxIter))
  coefficients <- out[-length(out)]
  names(coefficients) <- colnames(x)
  variance <- out[length(out)]
  names(variance) <- "var"

  # Fitting Random Effects
  matV <- matVFH(variance, samplingVar)
  if (is.null(x0Re)) {
    x0Re <- fitReStartingValues(y, x, coefficients, matV, psi)
  }
  re <- fitRe(
    y, x, coefficients, matV,
    x0Re, psi, addMaxIter(convCrit, maxIterRe)
  )

  # Iterations
  iterations <- storage$reformat(attr(out, "storage"))
  iterations <- c(iterations, re["iterations"])
  names(iterations) <- c("coefficients", "variance", "re")

  # Reformats
  re <- re$re
  reblup <- as.numeric(x %*% coefficients + matV$Z() %*% re)
  residuals <- y - reblup

  stripSelf(retList("fitrfh", public = c(
    "coefficients", "variance", "psi", "samplingVar", "y", "x", "iterations",
    "k", "tol", "K", "re", "reblup", "residuals"
  )))

}

fitCoefStartingValue <- function(y, x, matV) {
  xPrimeV <- crossprod(x, matV$VInv())
  solve(xPrimeV %*% x) %*% xPrimeV %*% y
}

fitReStartingValues <- function(y, x, beta, matV, psi) {
  # y: (numeric) response
  # x: (Matrix) Design-Matrix
  # beta: (numeric) estimated fixed effects
  # matV: (list) list of functions see e.g. matVFH
  U <- matU(matV$V())
  resids <- psi(U$sqrtInv() %*% (y - x %*% beta))
  as.numeric(
    tcrossprod(matV$Vu(), matV$Z()) %*% matV$VInv() %*% U$sqrt() %*% psi(resids)
  )
}

fitRe <- function(y, x, beta, matV, x0, psi, convCrit) {
  # y: (numeric) response
  # x: (Matrix) Design-Matrix
  # beta: (numeric) estimated fixed effects
  # matV: (list) list of functions see e.g. matVFH
  # psi: (function) influence function
  # convCrit: (function) convergence criterion

  fpFun <- fixedPointRobustRandomEffect(
    y, x, beta, matV, psi
  )

  re <- fixedPoint(addHistory(fpFun), x0, convCrit)
  iterRe <- attr(re, "history")

  list(re = as.numeric(re), iterations = cbind(iterRe, i = 1:NROW(iterRe)))

}

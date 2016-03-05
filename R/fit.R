#' @rdname rfh
#' @export
fitrfh <- function(y, x, samplingVar, ...) {
  # Non interactive fitting function for robust FH

  fixedPointParam <- function(parent = parent.frame()) {
    # To connect the returned function to the calling environment:
    enclosingEnv <- environment()
    parent.env(enclosingEnv) <- parent
    function(param) {
      # Defines one Iteration in algorithm:
      param <- as.numeric(param)
      beta <- param[-length(param)]
      sigma2 <- param[length(param)]

      fpBeta <- addHistory(fixedPointRobustBeta(
        y, x, matVFun(sigma2), psi = psi
      ))

      beta <- fixedPoint(fpBeta, beta, addMaxIter(convCrit, maxIter))

      fpSigma2 <- fixedPointRobustDelta(
        y, x, beta, matVFun, psi, K
      ) %>% addConstraintMin(0) %>% addHistory

      sigma2 <- fixedPoint(fpSigma2, sigma2, addMaxIter(convCrit, maxIter))

      list(beta, sigma2)
    }
  }

  matVFun <- . %>% matVFH(.samplingVar = samplingVar)
  out <- fitGenericModel(y, x, matVFun, fixedPointParam, ...)
  names(out$iterations) <- c("coefficients", "variance", "re")
  names(out$variance) <- "var"

  stripSelf(retList("fitrfh", public = c("samplingVar"), super = out))

}

fitGenericModel <- function(
  y, x, matVFun, fixedPointParam,
  k = 1.345, K = getK(k), psi = . %>% psiOne(k),
  x0Coef = NULL, x0Var = 1, x0Re = NULL,
  tol = 1e-6, maxIter = 100, maxIterRe = 100, convCrit = convCritRelative(tol)) {
  # Non interactive fitting function for robust FH Models

  # Fitting Model Parameter:
  if (is.null(x0Coef)) {
    x0Coef <- as.numeric(fitCoefStartingValue(y, x, matVFun(x0Var)))
  }
  out <- fixedPoint(
    addStorage(fixedPointParam()),
    c(x0Coef, x0Var),
    addMaxIter(convCrit, maxIter)
  )
  coefficients <- out[1:length(x0Coef)]
  names(coefficients) <- colnames(x)
  variance <- out[(length(x0Coef) + 1):length(out)]

  # Fitting Random Effects
  matV <- matVFun(variance)
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

  # Reformats
  re <- re$re
  reblup <- as.numeric(x %*% coefficients + matV$Z() %*% re)
  residuals <- y - reblup

  stripSelf(retList(public = c(
    "coefficients", "variance", "iterations",
    "tol", "maxIter", "maxIterRe",
    "k", "K", "psi",
    "y", "x", "re", "reblup", "residuals"
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

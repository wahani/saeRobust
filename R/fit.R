#' Fitting Precedures
#'
#' Several fitting procedures. Not intended for interactive use.
#'
#' @param y (numeric) response vector
#' @param x ([m|M]atrix) the design matrix
#' @param samplingVar (numeric) vector with sampling variances
#' @param ... arguments passed to \code{fitGenericModel}
#' @param matVFun (function) a function with one argument - the variance
#'   parameters - constructing something like \link{matVFH}
#' @param fixedPointParam (function) a function with one argument. The vector of
#'   model parameters. Returns a list of results of the next iteration in the
#'   overall algorithm.
#' @param k (numeric) tuning constant
#' @param K (numeric) scaling constant
#' @param psi (function) influence function
#' @param x0Coef (numeric) starting values for regression coefficients
#' @param x0Var (numeric) starting values for variance parameters
#' @param x0Re (numeric) starting values for random effects
#' @param tol (numeric) numerical toloerance to be used during optimisation
#' @param maxIter (integer) the maximum number of iterations
#' @param maxIterRe (integer) the maximum number of iterations for fitting the
#'   random effects
#' @param convCrit (function) a function defining the stopping rule
#'
#' @rdname fit
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
  names(out$variance) <- "variance"

  stripSelf(retList("fitrfh", public = c("samplingVar"), super = out))

}

#' @rdname fit
#' @export
fitrsfh <- function(y, x, samplingVar, W, x0Var = c(0.5, 1), ...) {
  # Non interactive fitting function for robust FH

  fixedPointParam <- function(parent = parent.frame()) {
    # To connect the returned function to the calling environment:
    enclosingEnv <- environment()
    parent.env(enclosingEnv) <- parent
    function(param) {
      # Defines one Iteration in algorithm:
      param <- as.numeric(param)
      beta <- param[-c(length(param)-1, length(param))]
      sigma2 <- param[length(param)]
      rho <- param[length(param) - 1]

      fpBeta <- addHistory(fixedPointRobustBeta(
        y, x, matVFun(c(rho, sigma2)), psi = psi
      ))

      beta <- fixedPoint(fpBeta, beta, addMaxIter(convCrit, maxIter))

      fpSigma2 <- fixedPointRobustDelta(
        y, x, beta, function(x) matVFun(c(rho, x), "sigma2"), psi, K
        ) %>% addConstraintMin(0) %>% addHistory

      sigma2 <- fixedPoint(fpSigma2, sigma2, addMaxIter(convCrit, maxIter))

      fpRho <- fixedPointRobustDelta(
        y, x, beta, function(x) matVFun(c(x, sigma2), "rho"), psi, K
      ) %>%
        addAverageDamp %>%
        addConstraintMin(-0.99) %>% addConstraintMax(0.99) %>%
        addHistory

      rho <- fixedPoint(fpRho, rho, addMaxIter(convCrit, maxIter))

      list(beta, rho, sigma2)
    }
  }

  matVFun <- function(x, ...) {
    matVSFH(.rho = x[1], .sigma2 = x[2], .W = W, .samplingVar = samplingVar, ...)
  }
  out <- fitGenericModel(y, x, matVFun, fixedPointParam, x0Var = x0Var, ...)
  names(out$iterations) <- c("coefficients", "correlation", "variance", "re")
  names(out$variance) <- c("correlation", "variance")

  stripSelf(retList(
    c("fitrsfh", "fitrfh"),
    public = c("samplingVar", "W"), super = out)
  )

}

#' @export
#' @rdname fit
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

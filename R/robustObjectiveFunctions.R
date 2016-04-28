#' Robust score function (ML) for beta
#'
#' Constructs a list of functions with \code{f} as the score and \code{f1} as
#' its derivative. Both are functions of the beta coefficients.
#'
#' @param y vector of response
#' @param x design matrix
#' @param matV (list of functions) see \link{matVFH}
#' @param psi influence function
#'
#' @rdname objectiveFunctions
#'
#' @export
scoreRobustBeta <- function(y, x, matV, psi) {
    # Helper functions
    resid <- function(beta) U$sqrtInv() %*% (y - x %*% beta)
    D <- function(beta) Diagonal(x = psi(resid(beta), deriv = TRUE))

    # Precalculations - they only have to be done once
    U <- matU(matV$V())
    memP0 <- crossprod(x, matV$VInv())
    memP1 <- memP0 %*% U$sqrt()

    f <- function(beta) memP1 %*% psi(resid(beta))
    f1 <- function(beta) - memP0 %*% D(beta) %*% x

    list(f = f, f1 = f1)
}

#' Fixed Point Functions
#'
#' This is an implementation of robustified fixed point functions to identify
#' model parameters in any mixed linear model: regression coefficients, variance
#' parameters, and random effects
#'
#' @param y vector of response
#' @param x design matrix
#' @param matV (list of functions) see \link{matVFH}
#' @param psi influence function
#' @param stepSize (numeric) size to be used in numeric derivative
#' @param lowerBound (numeric) a lower bound, such that \code{param - stepSize}
#'   cannot be outside of the parameter space
#'
#' @rdname fixedPointFunctions
#' @export
fixedPointRobustBeta <- function(y, x, matV, psi) {
    makeMatA <- matAConst(y, x, matV, psi)
    function(beta) {
        as.numeric(makeMatA(beta) %*% y)
    }
}

#' @export
#' @rdname fixedPointFunctions
robustObjectiveDelta <- function(y, x, beta, matVFun, psi, K, derivSelect) {
  # This is the squared estimation equation for a variance parameter. It can be
  # used when the fixed point for delta
  function(rho) {

    matV <- matVFun(rho)
    U <- matU(matV$V())
    psiResid <- psi(U$sqrtInv() %*% (y - x %*% beta))

    as.numeric(
      crossprod(psiResid, U$sqrt()) %*% matV$VInv() %*%
        matV$deriv[[derivSelect]]() %*%
        matV$VInv() %*% U$sqrt() %*% psiResid -
        matTrace(K * matV$VInv() %*% matV$deriv[[derivSelect]]())
    )
  }
}

#' @export
#' @rdname fixedPointFunctions
fixedPointNumericDelta <- function(y, x, beta, matVFun, psi, K, derivSelect, stepSize, lowerBound) {
  obDelta <- robustObjectiveDelta(y, x, beta, matVFun, psi, K, derivSelect)
  function(rho) {
    rho + obDelta(rho) / ((obDelta(max(lowerBound, rho - stepSize)) - obDelta(rho)) / stepSize)
  }
}

#' @param matVFun a function with one argument constructing something similar to
#'   \link{matVFH}: e.g. \code{. \%>\% matVFH(c(1, 1))}
#' @param K constant, see \link{getK}
#' @param beta beta coefficients to be used
#' @param derivSelect an indicator to select the derivative in matV. Position or
#'   name
#'
#' @rdname fixedPointFunctions
#' @export
fixedPointRobustDelta <- function(y, x, beta, matVFun, psi, K, derivSelect = 1) {
    # Precalculations - they only have to be done once
    mem1 <- (y - x %*% beta)

    function(param) {
        matV <- matVFun(param)
        U <- matU(matV$V())
        resid <- U$sqrtInv() %*% mem1
        psiResid <- psi(resid)
        c1 <- K / param * matTrace(matV$VInv() %*% matV$deriv[[derivSelect]]())
        c2 <- crossprod(psiResid, U$sqrt()) %*% matV$VInv() %*%
            matV$deriv[[derivSelect]]() %*% matV$VInv() %*% U$sqrt() %*% psiResid

        as.numeric(c2 / c1)
    }
}

#' @rdname fixedPointFunctions
#' @export
fixedPointRobustDelta2 <- function(y, x, beta, matVFun, psi, K, derivSelect) {
  # Precalculations - they only have to be done once
  mem1 <- (y - x %*% beta)

  function(param) {
    matV <- matVFun(param)
    U <- matU(matV$V())
    resid <- U$sqrtInv() %*% mem1
    psiResid <- psi(resid)

    C1tmp <- K * matV$VInv() %*% matV$deriv[[derivSelect[1]]]() %*% matV$ZVuZInv()
    C2tmp <- K * matV$VInv() %*% matV$deriv[[derivSelect[2]]]() %*% matV$ZVuZInv()

    C1 <- matrix(ncol = 2, c(
      matTrace(C1tmp %*% matV$ZVuBarZ[[derivSelect[1]]]()),
      matTrace(C1tmp %*% matV$ZVuBarZ[[derivSelect[2]]]()),
      matTrace(C2tmp %*% matV$ZVuBarZ[[derivSelect[1]]]()),
      matTrace(C2tmp %*% matV$ZVuBarZ[[derivSelect[2]]]())
    ))

    c2 <- lapply(matV$deriv[derivSelect], function(deriv) {
      as.numeric(
        crossprod(psiResid, U$sqrt()) %*% matV$VInv() %*%
          deriv() %*% matV$VInv() %*% U$sqrt() %*% psiResid
      )}) %>% unlist

    as.numeric(solve(C1) %*% c2)
  }
}

#' @rdname fixedPointFunctions
#' @export
fixedPointRobustRandomEffect <- function(y, x, beta, matV, psi) {

    makeMatB <- matBConst(y, x, beta, matV, psi)
    memResid <- y - x %*% beta
    function(u) as.numeric(as.matrix(makeMatB(u)) %*% memResid)
}

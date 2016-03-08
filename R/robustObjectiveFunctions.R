#' Robust score function (ML) for beta
#'
#' Constructs a list of functions with \code{f} as the score and \code{f1} as
#' its derivative. Both are functions of the beta coefficients.
#'
#' @param y vector of response
#' @param X design matrix
#' @param matV (list of functions) see \link{matVFH}
#' @param psi influence function
#'
#' @rdname objectiveFunctions
#'
#' @export
scoreRobustBeta <- function(y, X, matV, psi) {
    # Helper functions
    resid <- function(beta) U$sqrtInv() %*% (y - X %*% beta)
    D <- function(beta) Diagonal(x = psi(resid(beta), deriv = TRUE))

    # Precalculations - they only have to be done once
    U <- matU(matV$V())
    memP0 <- crossprod(X, matV$VInv())
    memP1 <- memP0 %*% U$sqrt()

    f <- function(beta) memP1 %*% psi(resid(beta))
    f1 <- function(beta) - memP0 %*% D(beta) %*% X

    list(f = f, f1 = f1)
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
fixedPointNumericDelta <- function(y, x, beta, matVFun, psi, K, derivSelect, stepSize) {
  obRho <- robustObjectiveDelta(y, x, beta, matVFun, psi, K, derivSelect)
  function(rho) {
    rho + obRho(rho) / ((obRho(rho - stepSize) - obRho(rho)) / stepSize)
  }
}

#' Fixed Point Functions
#'
#' This is an implementation of a robustified fixed point function to identify
#' beta coefficients in any mixed linear model.
#'
#' @param y vector of response
#' @param X design matrix
#' @param matV (list of functions) see \link{matVFH}
#' @param psi influence function
#'
#' @rdname fixedPointFunctions
#' @export
fixedPointRobustBeta <- function(y, X, matV, psi) {
    makeMatA <- matAConst(y, X, matV, psi)
    function(beta) {
        as.numeric(makeMatA(beta) %*% y)
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
fixedPointRobustDelta <- function(y, X, beta, matVFun, psi, K, derivSelect = 1) {
    # Precalculations - they only have to be done once
    mem1 <- (y - X %*% beta)

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
fixedPointRobustDelta2 <- function(y, X, beta, matVFun, psi, K, derivSelect) {
  # Precalculations - they only have to be done once
  mem1 <- (y - X %*% beta)

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
fixedPointRobustRandomEffect <- function(y, X, beta, matV, psi) {
    makeMatB <- matBConst(y, X, beta, matV, psi)
    memResid <- y - X %*% beta
    function(u) as.numeric(makeMatB(u) %*% memResid)
}

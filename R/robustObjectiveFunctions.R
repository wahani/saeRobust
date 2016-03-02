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
#'
#' @rdname fixedPointFunctions
#' @export
fixedPointRobustDelta <- function(y, X, beta, matVFun, psi, K) {
    # Precalculations - they only have to be done once
    mem1 <- (y - X %*% beta)

    function(param) {
        matV <- matVFun(param)
        U <- matU(matV$V())
        resid <- U$sqrtInv() %*% mem1
        psiResid <- psi(resid)
        c1 <- matTrace(K / param * matV$VInv() %*% matV$deriv[[1]]())
        c2 <- crossprod(psiResid, U$sqrt()) %*% matV$VInv() %*%
            matV$deriv[[1]]() %*% matV$VInv() %*% U$sqrt() %*% psiResid

        as.numeric(c2 / c1)
    }
}

#' @rdname fixedPointFunctions
#' @export
fixedPointRobustRandomEffect <- function(y, X, beta, matV, psi) {
    makeMatB <- matBConst(y, X, beta, matV, psi)
    memResid <- y - X %*% beta
    function(u) as.numeric(makeMatB(u) %*% memResid)
}

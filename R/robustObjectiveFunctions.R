#' Robust score function (ML) for beta
#'
#' Constructs a list of functions with \code{f} as the score and \code{f1} as
#' its derivative. Both are functions of the beta coefficients.
#'
#' @param y vector of response
#' @param X design matrix
#' @param V variance matrix
#' @param Vinv inverse of V (optional)
#' @param psi influence function
#' @param resid function to compute residuals (optional)
#'
#' @export
scoreRobustBeta <- function(y, X, V, Vinv = solve(V), psi, resid = NULL) {
    force(y); force(psi)

    # Helper functions
    if (is.null(resid)) resid <- function(beta) U$sqrtInv %*% (y - X %*% beta)
    D <- function(beta) Diagonal(x = psi(resid(beta), deriv = TRUE))

    # Precalculations - they only have to be done once
    U <- matU(V)
    memP0 <- crossprod(X, Vinv)
    memP1 <- memP0 %*% U$sqrt

    f <- function(beta) memP1 %*% psi(resid(beta))
    f1 <- function(beta) - memP0 %*% D(beta) %*% X

    list(f = f, f1 = f1)
}

#' Fixed Point Functions
#'
#' This is an implementation of a robustified fixed point function to identify
#' beta coefficients in any mixed linear model. The function is derived from the
#' pseudo linear form of robust mixed models in CCST (2015) and CCT (2011).
#'
#' @param y vector of response
#' @param X design matrix
#' @param V variance matrix
#' @param Vinv inverse of V (optional)
#' @param psi influence function
#' @param resid function to compute residuals (optional)
#'
#' @rdname fixedPointFunctions
#' @export
fixedPointRobustBeta <- function(y, X, V, Vinv = solve(V), psi, resid = NULL) {
    force(y)
    makeMatA <- matAConst(y, X, V, Vinv, psi)
    function(beta) {
        as.numeric(makeMatA(beta) %*% y)
    }
}

#' @param samplingVar the sampling variance of the direct estimator.
#' @param K constant, see \link{getK}. \code{length(beta) == ncols(X)}
#' @param beta beta coefficients to be used.
#'
#' @rdname fixedPointFunctions
#' @export
fixedPointRobustVarianceFH <- function(y, X, samplingVar, psi, K, beta) {
    force(samplingVar); force(psi); force(K)

    # Precalculations - they only have to be done once
    mem1 <- (y - X %*% beta)

    function(sigma2) {
        V <- matVFH(sigma2, samplingVar)
        U <- matU(V$V)
        resid <- U$sqrtInv %*% mem1
        psiResid <- psi(resid)

        a <- t(psiResid) %*% U$sqrt %*% V$vInv %*% V$vInv %*% U$sqrt %*% psiResid
        A <- matTrace(K * V$vInv %*% V$gInv)

        max(0, as.numeric(a / A))
    }
}

#' @param VuSqrtInv (Matrix) n times n matrix of the variance structure of the
#'   random effects part to the power of -1/2
#' @param VeSqrtInv (Matrix) n times n matrix of the variance structure of the
#'   model error part to the power of -1/2
#'
#' @rdname fixedPointFunctions
#' @export
fixedPointRobustRandomEffect <- function(y, X, beta, VuSqrtInv, VeSqrtInv, psi) {
    makeMatB <- matBConst(y, X, beta, VuSqrtInv, VeSqrtInv, psi)
    memResid <- y - X %*% beta
    function(u) as.numeric(makeMatB(u) %*% memResid)
}

#' Robust score function (ML) for beta
#'
#' Constructs a list of functions with \code{f} as the score and \code{f1} as its derivative. Both are functions of the beat coefficients.
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
    if(is.null(resid)) resid <- function(beta) sqrtUinv %*% (y - X %*% beta)
    D <- function(beta) diag(psi(resid(beta), deriv = TRUE))

    # Precalculations:
    # Sinha & Rao (2009): page 386 for the def of U
    U <- diag(diag(V))
    sqrtU <- sqrt(U)
    sqrtUinv <- diag(1 / diag(sqrtU))

    memP0 <- crossprod(X, Vinv)
    memP1 <- memP0 %*% sqrtU

    f <- function(beta) memP1 %*% psi(resid(beta))
    f1 <- function(beta) - memP0 %*% D(beta) %*% X

    list(f = f, f1 = f1)
}

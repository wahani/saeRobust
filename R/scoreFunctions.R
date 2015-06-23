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
    if(is.null(resid)) resid <- function(beta) U$sqrtInv %*% (y - X %*% beta)
    D <- function(beta) Diagonal(x = psi(resid(beta), deriv = TRUE))

    # Precalculations - they only have to be done once
    U <- matU(V)
    memP0 <- crossprod(X, Vinv)
    memP1 <- memP0 %*% U$sqrt

    f <- function(beta) memP1 %*% psi(resid(beta))
    f1 <- function(beta) - memP0 %*% D(beta) %*% X

    list(f = f, f1 = f1)
}

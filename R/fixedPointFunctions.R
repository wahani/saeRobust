#' Fixed Point Function for Robust Beta
#'
#' This is an implementation of a robustified fixed point function to identify beta coefficients in any mixed linear model. The function is derived from the pseudo linear form of robust mixed models in CCST (2015) and CCT (2011).
#'
#' @param y vector of response
#' @param X design matrix
#' @param V variance matrix
#' @param Vinv inverse of V (optional)
#' @param psi influence function
#' @param resid function to compute residuals (optional)
#'
#' @export
fixedPointRobustBeta <- function(y, X, V, Vinv = solve(V), psi, resid = NULL) {
    force(y); force(psi)

    # Helper functions
    if(is.null(resid)) resid <- function(beta) U$sqrtInv %*% (y - X %*% beta)
    wj <- function(beta) {
        resids <- resid(beta)
        psi(resids) / resids
    }

    # Precalculations - they only have to be done once
    U <- matU(V)
    memP0 <- crossprod(X, Vinv)
    memP1 <- memP0 %*% U$sqrt
    memP2 <- U$sqrtInv %*% X

    function(beta) {
        W <- Diagonal(x = wj(beta))
        as.numeric(solve(memP1 %*% W %*% memP2) %*% memP1 %*% W %*% U$sqrtInv %*% y)
    }
}

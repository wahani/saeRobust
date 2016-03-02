#' Matrix constructor functions
#'
#' @param .V (Matrix) variance matrix
#'
#' @details
#'
#' \code{matU} computes U. U is the matrix containing only the diagonal
#'   elements of V. This function returns a list of functions which can be
#'   called to compute specific transformations of U.
#'
#' @rdname varianceMatrices
#' @export
matU <- function(.V) {

    .diag <- function(x) Diagonal(x = x)

    U <- getter(diag(.V), .diag)
    sqrt <- getter(base::sqrt(diag(.V)), .diag)
    sqrtInv <- getter(1 / base::sqrt(diag(.V)), .diag)
    retList()

}

#' @param x ([m|M]atrix) a matrix
#'
#' @details \code{matTrace} computes the trace of a matrix.
#'
#' @rdname varianceMatrices
#' @export
matTrace <- function(x) {
    sum(diag(x))
}

#' @param .sigma2 (numeric) a scalar. The variance parameter of the random
#'   effects.
#' @param .samplingVar (numeric) the vector of sampling variances of the direct
#'   estimator.
#'
#' @details \code{matVFH} constructs variance-covariance matrix for a FH-model.
#'   Returns a list of functions to compute various transformations.
#'
#' @rdname varianceMatrices
#' @export
matVFH <- function(.sigma2, .samplingVar) {

    .diag <- function(x) Diagonal(x = x)

    Vu <- getter(rep(.sigma2, length(.samplingVar)), .diag)
    VuInv <- getter(1 / rep(.sigma2, length(.samplingVar)), .diag)
    Ve <- getter(.samplingVar, .diag)
    VeInv <- getter(1 / .samplingVar, .diag)
    V <- getter(Vu() + Ve())
    VInv <- getter(solve(V()))
    Z <- getter(Diagonal(length(.samplingVar)))

    deriv <- list(
        getter(Diagonal(length(.samplingVar)))
    )

    retList()

}

#' @param y (numeric) response
#' @param X (Matrix) design matrix
#' @param beta (numeric) vector of regression coefficients
#' @param u (numeric) vector of random effects
#' @param matV (list of functions) see \code{matVFH} for an example
#' @param psi (function) the influence function
#'
#' @details \code{matB} computes the matrix B which is used to compute the
#'   weights in the pseudo linearised representation of the REBLUP.
#'
#' @rdname varianceMatrices
#' @export
matB <- function(y, X, beta, u, matV, psi) {
    matBConst(y, X, beta, matV, psi)(u)
}

#' @details \code{matBConst} returns a function with one argument, u, to compute
#'   the matrix B. This function is used internally to compute B in the fixed
#'   point algorithm (see \link{fixedPointRobustRandomEffect}).
#'
#' @rdname varianceMatrices
#' @export
matBConst <- function(y, X, beta, matV, psi) {

    Ue <- matU(matV$Ve())
    Uu <- matU(matV$Vu())

    # Helper functions
    resid <- function(u) as.numeric(memResid - u)

    w2 <- function(u) {
        resids <- resid(u) * diag(Ue$sqrtInv())
        psi(resids) / resids
    }

    w3 <- function(u) {
        resids <- u * diag(Uu$sqrtInv())
        psi(resids) / resids
    }

    # Precalculations - they only have to be done once
    memXB <- X %*% beta
    memResid <- y - memXB
    memZVeU <- crossprod(matV$Z(), matV$VeInv())

    function(u) {
        W2 <- Diagonal(x = w2(u))
        W3 <- Diagonal(x = w3(u))
        memPart1 <- memZVeU %*% W2
        memPart2 <- matV$VuInv() %*% W3
        solve(memPart1 %*% matV$Z() + memPart2) %*% memPart1
    }
}

#' @details \code{matA} computes the matrix A which is used to compute the
#'   weights in the pseudo linearized representation of the REBLUP.
#'
#' @rdname varianceMatrices
#' @export
matA <- function(y, X, beta, matV, psi) {
    matAConst(y, X, matV, psi)(beta)
}

#' @details \code{matAConst} returns a function with one argument, beta, to
#'   compute the matrix A. This function is used internally to compute A in the
#'   fixed point algorithm for beta (see \link{fixedPointRobustBeta}).
#'
#' @rdname varianceMatrices
#' @export
matAConst <- function(y, X, matV, psi) {
    force(y); force(psi)

    # Helper functions
    resid <- function(beta) as.numeric(U$sqrtInv() %*% (y - X %*% beta))

    w <- function(beta) {
        resids <- resid(beta)
        psi(resids) / resids
    }

    # Precalculations - they only have to be done once
    U <- matU(matV$V())
    memP0 <- crossprod(X, matV$VInv())

    function(beta) {
        W1 <- Diagonal(x = w(beta))
        solve(memP0 %*% W1 %*% X) %*% memP0 %*% W1
    }
}

#' @details \code{matW} returns a matrix containing the weights as they are
#'   defined for the pseudo linear form, such that \code{matW \%*\% y} is the
#'   REBLUP.
#'
#' @rdname varianceMatrices
#' @export
matW <- function(y, X, beta, u, matV, psi) {
    A <- matA(y, X, beta, matV, psi)
    B <- matB(y, X, beta, u, matV, psi)
    XA <- X %*% A
    XA + B %*% (Diagonal(length(y)) - XA)
}

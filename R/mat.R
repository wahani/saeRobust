#' Helper functions operating on matrices and matrix constructors
#'
#' @param V variance matrix
#'
#' @details \code{matU} computes U. U is the matrix containing only the diagonal elements of V. See Sinha & Rao (2009): page 386 for the def of U.
#'
#' @rdname varianceMatrices
#' @export
matU <- function(V) {
    U <- Diagonal(x = diag(V))
    sqrt <- sqrt(U)
    sqrtInv <- solve(sqrt)
    retList()
}

#' @param sigma2 a scalar. The variance parameter of the random effects.
#' @param samplingVar the vector of sampling variances of the direct estimator.
#'
#' @details \code{matVFH} constructs variance-covariance matrix for a FH-model. Returns a list with the matrix and its inverse.
#'
#' @rdname varianceMatrices
#' @export
matVFH <- function(sigma2, samplingVar) {
    G <- Diagonal(length(samplingVar), sigma2)
    gInv <- solve(G)
    R <- Diagonal(x = samplingVar)
    rInv <- solve(R)
    V <- G + R
    vInv <- solve(V)
    retList()
}

#' @param x a matrix
#'
#' @details \code{matTrace} computes the trace of a matrix.
#'
#' @rdname varianceMatrices
#' @export
matTrace <- function(x) {
    sum(diag(x))
}

#' Computes the U matrix
#'
#' Sinha & Rao (2009): page 386 for the def of U
#'
#' @param V variance matrix
#'
#' @export
matU <- function(V) {
    U <- Diagonal(x = diag(V))
    sqrtU <- sqrt(U)
    sqrtUinv <- solve(sqrtU)

    list(U = U, sqrt = sqrtU, sqrtInv = sqrtUinv)
}

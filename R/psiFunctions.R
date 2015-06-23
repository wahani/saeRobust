#' psiOne
#'
#' @description influence-function
#'
#' @param u standardized residuals
#' @param k tuning constant
#' @param deriv if \code{TRUE} returns the derivative
#'
#' @rdname psiFunctions
#' @export
psiOne <- function(u, k = 1.345, deriv = FALSE){
    u <- as.numeric(u) # in the case that u comes from the Matrix package

    if(deriv) return(as.numeric(abs(u) <= k))

    sm <- median(abs(u)) / 0.6745
    w <- pmin(1, k / abs(u / sm))
    w * u
}

#' @rdname psiFunctions
#' @export
getK <- function(k) {
    2 * pnorm(k) - 1 - 2 * k * dnorm(k) + 2 * (k^2) * (1-pnorm(k))
}

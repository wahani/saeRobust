#' psiOne
#'
#' @description influence-function
#'
#' @param u standardized residuals
#' @param k tuning constant
#' @param deriv if \code{TRUE} returns the derivative
#'
#' @export
psiOne <- function(u, k = 1.345, deriv = FALSE){

    if(deriv) return(as.numeric(abs(u) <= k))

    sm <- median(abs(u)) / 0.6745
    w <- pmin(1, k / abs(u / sm))
    w * u
}

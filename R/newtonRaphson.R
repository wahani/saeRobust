#' Newton-Raphson Algorithm
#'
#' @description Newton-Raphson finds zeroes of a function. The user can supply the function and its first derivative. Note that the Newton Raphson Algorithm is a special case of a fixed point algorithm thus it is implemented using \code{\link{fixedPoint}} and is only a convenience.
#'
#' @param funList the function to be evaluated in the algorithm
#' @param ... arguments passed to \code{\link{fixedPoint}}
#'
#' @export
#' @rdname newtonRaphson
#' @examples
#' \dontrun{
#' vignette("fixedPoint", "saeRobustTools")
#' }
newtonRaphson <- function(funList, ...) {
    force(funList)
    fun <- function(x) as.numeric(x - solve(funList$f1(x)) %*% funList$f(x))
    fixedPoint(fun, ...)
}

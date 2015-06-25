#' Fixed Point Algorithm Infrastructure
#'
#' @description A fixed-point function supplied by the user is iteratively evaluated. \code{averageDamp} can be used to add average damping to the function - this may have a positive effect on the speed of convergence.
#'
#' @param fun the function to be evaluated in the algorithm
#' @param x0 starting value
#' @param convCrit a function returning a logical scalar. Is called with two arguments; the first is the value from iteration n; the second is the value from iteration n-1
#' @param tolerance a numeric value > 0
#'
#' @export
#' @rdname fixedPoint
#' @examples
#' \dontrun{
#' vignette("fixedPoint", "saeRobustTools")
#' }
fixedPoint <- function(fun, x0, convCrit) {
    x1 <- NULL
    repeat {
        x1 <- fun(x0)
        if(convCrit(x1, x0)) {
            break
        } else {
            x0 <- x1
        }
    }
    x1
}

#' @export
#' @rdname fixedPoint
averageDamp <- function(fun) {
    force(fun)
    function(x) (x + fun(x)) / 2
}

#' @export
#' @rdname fixedPoint
convCritAbsolute <- function(tolerance = 1e-5) {
    assert_that(tolerance > 0)
    function(xn1, xn0) all(abs(xn0 - xn1) < tolerance)
}

#' @export
#' @rdname fixedPoint
convCritRelative <- function(tolerance = 1e-5) {
    assert_that(tolerance > 0)
    function(xn1, xn0) all((abs(xn0 - xn1) / max(1, abs(xn0))) < tolerance)
}


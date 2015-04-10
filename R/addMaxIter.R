#' Higher Order Function to Add a Maximum Iteration
#'
#' This function can be used to modify convergence criterion functions.
#'
#' @param fun
#' @param maxIter
#'
#' @export
addMaxIter <- function(fun, maxIter) {
    force(fun); force(maxIter)
    count <- 0
    function(...) {
        count <<- count + 1
        if (count > maxIter) TRUE else fun(...)
    }
}

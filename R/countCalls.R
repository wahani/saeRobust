#' Higher Order Function to Count the Number of Calls
#'
#' This function can be used to count the number of calls of a function. See \code{\link{fixedPoint}} for how it can be used.
#'
#' @param fun function
#'
#' @export
#' @rdname countCalls
countCalls <- function(fun) {
    force(fun)
    count <- 0
    function(...) {
        count <<- count + 1
        DataWithCount(fun(...), count = count)
    }
}

#' @exportClass DataWithCount
#' @rdname countCalls
DataWithCount <- setClass("DataWithCount", slots = c(.Data = "ANY", count = "numeric"))



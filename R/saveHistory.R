#' Higher Order Function to Save the History of Results
#'
#' This function can be used to save a history of results of a function. The history is stored as a matrix, so this works best if the return value of \code{fun} is numeric. See \code{\link{fixedPoint}} for how it can be used.
#'
#' @param fun function
#'
#' @export
#' @rdname saveHistory
saveHistory <- function(fun) {
    force(fun)
    history <- NULL
    function(...) {
        res <- fun(...)
        history <<- rbind(history, res)
        DataWithHistory(res, history = history)
    }
}

#' @exportClass DataWithHistory
#' @rdname saveHistory
DataWithHistory <- setClass("DataWithHistory", slots = c(.Data = "ANY", history = "matrix"))






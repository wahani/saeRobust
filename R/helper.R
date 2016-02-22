addAttr <- function(to, what, name) {
    attr(to, name) <- what
    to
}

getter <- function(.x, .wrapper = identity) {
    memoise::memoise(function() {
        .wrapper(.x)
    })
}

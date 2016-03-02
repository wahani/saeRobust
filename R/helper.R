addAttr <- function(to, what, name) {
    attr(to, name) <- what
    to
}

getter <- function(.x, .wrapper = identity) {
    memoise::memoise(function() {
        .wrapper(.x)
    })
}

storage <- module({
    reformat <- function(storage) {
        .iter <- function(ind) {
            map(storage, ~ attr(.[[ind]], "history")) %>%
            { map(seq_along(.), i ~ cbind(.[[i]], i)) } %>% # add iteration
                do.call(what = rbind)
        }

        map(seq_along(storage[[1]]), .iter)
    }
})

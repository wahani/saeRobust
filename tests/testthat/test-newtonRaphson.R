context("Newton Raphson")

test_that("Fixed Point Framework", {

    sqrtNr <- function(.p) {
        force(.p)
        f <- function(x) x^2 - .p
        f1 <- function(x) 2 * x
        as.list(environment())
    }

    nr <- function(...) newtonRaphson(..., convCrit = function(xn1, xn0) abs(xn0 - xn1) < 0.001)
    expect_equal(nr(sqrtNr(2), 2), sqrt(2), tolerance = 1e-3)
})

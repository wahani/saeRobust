context("Fixed Point")

test_that("Fixed Point Framework", {

    expect_equal <- testthat::expect_equal

    sqrtFp <- function(p) {
        force(p)
        function(x) p / x
    }

    fpsqrt2Damped <- averageDamp(sqrtFp(2))
    fp <- function(...) fixedPoint(..., convCrit = function(xn1, xn0) abs(xn0 - xn1) < 0.001)

    expect_equal(sqrtFp(2)(sqrt(2)), sqrt(2))
    expect_equal(fpsqrt2Damped(2), 1.5)
    funWithCounter <- countCalls(function() 1)
    funWithCounter(); funWithCounter()
    expect_equal(slot(funWithCounter(), "count"), 3)

    res <- fp(countCalls(fpsqrt2Damped), 2)
    expect_equal(slot(res, "count"), 4)
    expect_equal(res@.Data, sqrt(2))

    res <- fp(saveHistory(fpsqrt2Damped), 2)
    expect_equal(nrow(slot(res, "history")), 4)
    expect_equal(res@.Data, sqrt(2))

    res <- fp(addConstraintMin(fpsqrt2Damped, 1.5), 2)
    expect_equal(res, 1.5)

    res <- fp(addConstraintMax(fpsqrt2Damped, 1.2), 1)
    expect_equal(res, 1.2)

})

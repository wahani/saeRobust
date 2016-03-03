context("robust ee")
test_that("estimation equations are computed correctly", {

    set.seed(1)

    X <- cbind(1, 1:10)
    u <- rnorm(10)
    sigu <- var(u)
    y <- X %*% c(1, 1) + u

    matV <- matVFH(sigu, rep(0.01, 10))
    coefs <- .lm.fit(X, y)$coefficients
    uest <- .lm.fit(X, y)$residuals

    ree <- robEstEqu(y, X, coefs, uest, matV, identity, 1)

    testthat::expect_equal(ree$beta(), c(0, 0))

    # How to get the other vals to zero without fitting the model is beyond me:
    testthat::expect_is(ree$delta(), "numeric")
    testthat::expect_equal(length(ree$delta()), 1)

    testthat::expect_is(ree$re(), "numeric")
    testthat::expect_equal(length(ree$re()), 10)

})


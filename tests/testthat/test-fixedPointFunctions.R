context("Fixed-Point Functions")
test_that("Fixed Point for Robust Beta is correct", {
    # Test Data
    set.seed(1)
    y <- 2 * 1:5 + rnorm(100)
    X <- cbind(1, rep(1:5, 20))
    V <- diag(1, 100)
    Vinv <- solve(V)

    convCrit <- function(xn1, xn0) all(abs(xn0 - xn1) < 0.001)

    # Equivalence to Newton-Raphson
    scoreFuns <- scoreRobustBeta(y, X, V, Vinv, psiOne)
    fpFun <- fixedPointRobustBeta(y, X, V, Vinv, psiOne)

    solutionNR <- newtonRaphson(scoreFuns, x0 = c(0, 2), convCrit = convCrit)
    solutionFP <- fixedPoint(fpFun, c(0, 2), addMaxIter(convCrit, 2))

    expect_equal(solutionNR, solutionFP, tolerance = 1e-3)

    # equivalence to OLS
    fpFun <- fixedPointRobustBeta(y, X, V, Vinv, Curry(psiOne, k = 100))
    solutionFP <- fixedPoint(fpFun, c(0, 2), addMaxIter(convCrit, 2))
    expect_equal(solutionFP, as.numeric(solve(crossprod(X), crossprod(X, y))))
})


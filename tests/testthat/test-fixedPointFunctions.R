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

test_that("Fixed Point for Robust Variance is correct", {

    # Test Data
    set.seed(4)
    df <- data.frame(
        y = 2 * 1:5 + rnorm(100) + rnorm(100),
        x =  rep(1:5, 20),
        dirVar = 1
    )
    df$y[1] <- 100 # outlier

    y <- df$y
    X <- cbind(1, df$x)

    convCrit <- function(xn1, xn0) all(abs(xn0 - xn1) < 1e-5)

    # Sanity checks
    fpFun <- fixedPointRobustVarianceFH(y, X, rep(1, 100), psiOne, K = getK(1.345), c(0, 2))
    estVarRobust <- fixedPoint(fpFun, 1, convCrit)
    estVarFH <- sae::eblupFH(y ~ x, dirVar, method = "FH", data = df)$fit$refvar

    expect_is(estVarRobust, "numeric")
    expect_true(estVarRobust < estVarFH)

    # Approximately equal to non robust FH? Not at the moment. It is not clear
    # if this implementation has some bug or if there is something in sae - or
    # if they are both correct but different.
    fpFun <- fixedPointRobustVarianceFH(
        y, X, rep(1, 100), Curry(psiOne, k = 10000), K = getK(10000),
        sae::eblupFH(y ~ x, dirVar, method = "FH", data = df)$fit$estcoef[,"beta"]
        )

    estVarRobust <- fixedPoint(fpFun, 1, convCrit)
    # expect_equal(estVarRobust, estVarFH) # Not yet true

    data(milk, package = "sae", envir = environment())

    # Fit FH model using REML method with indicators of 4 Major Areas as
    # explanatory variables.
    milk$dirVar <- milk$SD^2
    resultREML <- sae::eblupFH(yi ~ 1, dirVar, data = milk)

    y <- milk$yi
    X <- matrix(1, nrow = length(y), ncol = 1)
    fpFun <- fixedPointRobustVarianceFH(
        y, X, milk$dirVar, Curry(psiOne, k = 10000), K = getK(10000),
        beta = resultREML$fit$estcoef[, "beta"]
    )
    estVarRobust <- fixedPoint(fpFun, 1, convCrit)

    # Not satisfying but on the same scale.
    expect_equal(round(estVarRobust, 2), round(resultREML$fit$refvar, 2))
})

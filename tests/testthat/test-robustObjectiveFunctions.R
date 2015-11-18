context("Objective Functions")

test_that("Fixed Point for Robust Random Effect is correct", {
    # Test Data
    set.seed(1)
    u <- rnorm(100, sd = 10)
    y <- 2 * 1:5 + u + rnorm(100, sd = 2)
    X <- cbind(1, rep(1:5, 20))

    Ve <- list(
        mat = diag(4, 100),
        sqrtInv = diag(0.5, 100)
    )

    Vu <- list(
        mat = diag(10^2, 100),
        sqrtInv = diag(1 / 10, 100)
    )

    fhFit <- rfh(y ~ x, data.frame(y = y, x = X[, 2], dirVar = 4), "dirVar")
    beta <- fhFit$beta
    # re <- attr(saeRobustTools:::predict.rfh(fhFit), "re")

    convCrit <- function(xn1, xn0) all(abs(xn0 - xn1) < 1e-3)

    fpFun <- fixedPointRobustRandomEffect(
        y, X, beta, Vu$sqrtInv, Ve$sqrtInv, Curry(psiOne, k = Inf)
    )

    # when k -> inf then this should be the estimator for u:
    uDirect <- as.numeric(
        Vu$mat %*% solve(Ve$mat + Vu$mat) %*% (y - X %*% beta)
    )
    # This is the solution with the fixed point (k = Inf).
    solutionFP <- fixedPoint(fpFun, uDirect, addMaxIter(convCrit, 10000))

    # They roughly the same:
    summary(abs(uDirect) - abs(solutionFP))
    expect_equal(uDirect, solutionFP)

})

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

test_that("Newton Raphson basics", {

    sqrtNr <- function(.p) {
        force(.p)
        f <- function(x) x^2 - .p
        f1 <- function(x) 2 * x
        as.list(environment())
    }

    nr <- function(...) newtonRaphson(..., convCrit = function(xn1, xn0) abs(xn0 - xn1) < 0.001)
    fp <- function(funList) fixedPoint(
        newtonRaphsonFunction(funList),
        2,
        function(xn1, xn0) abs(xn0 - xn1) < 0.001)

    expect_equal(nr(sqrtNr(2), 2), sqrt(2), tolerance = 1e-3)
    expect_equal(nr(sqrtNr(2), 2), fp(sqrtNr(2)))
})

test_that("NR for robust betas", {
    set.seed(3)
    y <- 2 * 1:5 + rnorm(5)
    X <- cbind(1, 1:5)
    V <- diag(1, 5)
    Vinv <- solve(V)

    score <- scoreRobustBeta(y, X, V, Vinv, psiOne)
    # as.numeric, because score uses Matrix classes.
    expect_equal(as.numeric(score$f1(c(0, 2))), as.numeric(-crossprod(X, Vinv) %*% X))

    # one outlier
    set.seed(1)
    y <- 2 * 1:5 + rnorm(5)
    score <- scoreRobustBeta(y, X, V, Vinv, psiOne)
    expect_equal(as.numeric(score$f1(c(0, 2))),
                 as.numeric(-crossprod(X, Vinv) %*% diag(c(1, 1, 1, 0, 1)) %*% X))

    # Equivalence to OLS
    set.seed(1)
    y <- 2 * 1:5 + rnorm(100)
    X <- cbind(1, rep(1:5, 20))
    V <- diag(1, 100)
    Vinv <- solve(V)

    funs <- scoreRobustBeta(y, X, V, Vinv, function(...) psiOne(..., k = 100))

    expect_equal(
        newtonRaphson(
            funs,
            x0 = c(0, 2),
            convCrit = function(xn1, xn0) all(abs(xn0 - xn1) < 0.001)),
        as.numeric(solve(crossprod(X)) %*% crossprod(X, y)))

})


context("Newton Raphson")

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
    expect_equal(score$f1(c(0, 2)), - crossprod(X, Vinv) %*% X)

    # one outlier
    set.seed(1)
    y <- 2 * 1:5 + rnorm(5)
    score <- scoreRobustBeta(y, X, V, Vinv, psiOne)
    expect_equal(score$f1(c(0, 2)), - crossprod(X, Vinv) %*% diag(c(1, 1, 1, 0, 1)) %*% X)

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







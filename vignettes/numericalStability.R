## ------------------------------------------------------------------------
library("Matrix")
library("magrittr")

# Test X
Z1 <- Matrix(rep(999, 4), ncol = 2)
X <- testMatX(Z1, Z1)
X
isSymmetric(X)
det(X)
solve(X)

# Test Y
y0 <- testResponse0(X)
y <- as.numeric(testResponse(y0, k = 4))
dirVar <- rep(var(as.numeric((abs(y - y0)))), length(y))

# fh <- sae::eblupFH(
#     y ~ X2 + X3 + X4 + X5 + X6 - 1,
#     dirVar,
#     data = data.frame(y = as.numeric(y), as.matrix(X), dirVar)
# )
# 
# fit <- fitrfh(y, X[, -1], dirVar, c(lm.fit(as.matrix(X[, -1]), y)$coefficients, 1), k = 1e6)
# sum((predict(fit)$REBLUP - y) / y)
# sum((fh$eblup - y) / y)



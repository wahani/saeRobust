## ------------------------------------------------------------------------
library("Matrix")
library("magrittr")

# Test X
Z1 <- Matrix(rep(998, 4), ncol = 2)
X <- testMatX(Z1, Z1)
X
isSymmetric(X)
det(X)
solve(X)

# Test Y
y0 <- testResponse0(X)
y <- testResponse(y0, k = 1)
dirVar <- var(as.numeric((abs(y - y0))))

lmModel <- lm(
    y ~ X2 + X3 + X4 + X5 + X6 - 1,
    data.frame(y = as.numeric(y), as.matrix(X))
)

tmp <- summary(lmModel)
tmp$residuals

fpBeta <- fixedPointRobustBeta(y, X, matVFH(0, rep(1, length(y))), psi = . %>% psiOne)
solutionFP <- fixedPoint(fpBeta, rep(1, ncol(X)), addMaxIter(convCritRelative(), 100))
y0Pred <- solutionFP %*% X %>% as.numeric
y0Pred - y

scoreFuns <- scoreRobustBeta(y, X, matVFH(0, rep(1, length(y))), psi = psiOne)
solutionNR <- newtonRaphson(
    scoreFuns, 
    x0 = rep(1, ncol(X)), 
    convCrit = addMaxIter(convCritRelative(), 100)
)

y0Pred <- solutionNR %*% X %>% as.numeric
y0Pred - y

solutionRFH <- saeRobustTools:::fitrfh(y, X[, -1], rep(dirVar, nrow(y)))

fhModel <- sae::eblupFH(
    y ~ X2 + X3 + X4 + X5 + X6 - 1,
    vardir,
    data = data.frame(y = as.numeric(y), as.matrix(X), vardir = dirVar)
)
 
sum((y - fhModel$eblup)^2) - sum((tmp$residuals)^2)




rfhfit <- function(y, X, samplingVar, theta0 = c(rep(1, ncol(X)), 1), convCrit = convCritAbsolute()) {
    # Non interactive fitting function for robust FH
    # y: (numeric) response
    # X: ((M|m)atrix) design matrix
    # samplingVar: (numeric) the sampling variances
    # theta0: (numeric) starting values
    # convCrit: (function) conversion criterion

    oneIter <- function(param) {
        beta <- param[-length(param)]
        sigma2 <- param[length(param)]

        fpBeta <- fixedPointRobustBeta(y, X, matVFH(sigma2, samplingVar)$V, psi = psiOne)
        beta <- fixedPoint(fpBeta, beta, convCrit)

        fpSigma2 <- fixedPointRobustVarianceFH(y, X, samplingVar, psiOne, getK(1.345), beta = beta)
        sigma2 <- fixedPoint(averageDamp(fpSigma2), sigma2, convCrit)

        as.numeric(c(beta, sigma2))
    }

    out <- fixedPoint(oneIter, theta0, convCrit)

    beta <- out[-length(out)]
    variance <- out[length(out)]
    retList("rfh", c("beta", "variance"))

}


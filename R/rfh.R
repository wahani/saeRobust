#' Robust Fay Herriot
#'
#' Implementations of robust FH models
#'
#' @param y the response
#' @param X the design matrix
#' @param samplingVar the vector of sampling variances
#' @param theta0 starting values for
#'
#' @rdname rfh
#' @export
rfh <- function(y, X, samplingVar, theta0 = c(rep(1, ncol(X)), 1)) {

    oneIter <- function(param) {
        # Defines one iteration of the overall algorithm
        beta <- param[-length(param)]
        sigma2 <- param[length(param)]

        fpBeta <- fixedPointRobustBeta(y, X, matVFH(sigma2, samplingVar)$V, psi = psiOne)
        beta <- fixedPoint(fpBeta, beta, convCritAbsolute())
        # beta <- fpBeta(beta)

        fpSigma2 <- fixedPointRobustVarianceFH(y, X, samplingVar, psiOne, getK(1.345), beta = beta)
        sigma2 <- fixedPoint(averageDamp(fpSigma2), sigma2, convCritAbsolute())
        # sigma2 <- averageDamp(fpSigma2)(sigma2)

        as.numeric(c(beta, sigma2))

    }

    out <- fixedPoint(printSteps((oneIter)), theta0, convCritAbsolute())

    beta <- out[-length(out)]
    variance <- out[length(out)]
    retList("rfh", c("beta", "variance"))

}

printSteps <- function(fun) {
    force(fun)
    function(...) {
        print(...)
        fun(...)
    }
}

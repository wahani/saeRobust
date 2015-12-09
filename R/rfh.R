#' Robust Fay Herriot Model
#'
#' @param formula (formula)
#' @param data (data.frame)
#' @param samplingVar (character)
#' @param ... arguments passed to methods
#'
#' @rdname rfh
#'
#' @export
rfh(formula, data, samplingVar, ...) %g% standardGeneric("rfh")

#' @rdname rfh
#' @export
rfh(formula ~ formula, data ~ data.frame, samplingVar ~ character, ...) %m% {
    call <- match.call()
    xy <- makeXY(formula, data)
    samplingVar <- data[[samplingVar]]

    retList(
        public = c("call", "xy", "samplingVar"),
        super = rfhfit(xy$y, xy$x, samplingVar)
    )

}

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

        fpBeta <- fixedPointRobustBeta(
            y, X, matVFH(sigma2, samplingVar)$V, psi = psiOne
        )
        beta <- fixedPoint(fpBeta, beta, convCrit)

        fpSigma2 <- fixedPointRobustVarianceFH(
            y, X, samplingVar, psiOne, getK(1.345), beta = beta
        )
        sigma2 <- fixedPoint(averageDamp(fpSigma2), sigma2, convCrit)

        as.numeric(c(beta, sigma2))
    }

    out <- fixedPoint(oneIter, theta0, convCrit)

    beta <- out[-length(out)]
    variance <- out[length(out)]
    retList("rfh", c("beta", "variance")) %>% stripSelf

}

#' @param object (rfh) an object of class rfh
#' @rdname rfh
#' @export
predict.rfh <- function(object, ...) {
    vC <- matVFH(object$variance, object$samplingVar)
    re <- reFitCCST(
        object$xy$y, object$xy$x, object$beta,
        vC$G, sqrt(vC$gInv),
        vC$R, sqrt(vC$rInv)
    )

    out <- as.numeric(object$xy$x %*% object$beta + re)
    attr(out, "re") <- as.numeric(re)
    out

}

reFitCCST <- function(y, X, beta, Vu, VuSqrtInv, Ve, VeSqrtInv,
                      psi = Curry(psiOne, k = 1.345),
                      convCrit = convCritAbsolute()) {
    # y: (numeric) response
    # X: (Matrix) Design-Matrix
    # beta: (numeric) estimated fixed effects
    # Vu: (Matrix) estimated variance-covariance of the random effect part
    # VuSqrtInv: (Matrix)
    # Ve: (Matrix) given sampling error Matrix
    # VeSqrtInv (Matrix)
    # psi (function) influence function

    # Non-robust random effects as starting values:
    startingValues <- as.numeric(Vu %*% solve(Ve + Vu) %*% (y - X %*% beta))

    fpFun <- fixedPointRobustRandomEffect(
        y, X, beta, VuSqrtInv, VeSqrtInv, psi
    )

    fixedPoint(fpFun, startingValues, convCrit)
}

# this is my old version. It uses a Newton-Raphson algorithm which has a
# potential bug. I need to reimpement and test that algorithm within this
# package
#
# predict.rfh <- function(object, ...) {
#     # This interface should be replaced in time.
#     interfaceList <- list(
#         reVar = object$variance,
#         vardir = object$samplingVar,
#         y = object$xy$y,
#         X = as.matrix(object$xy$x),
#         beta = object$beta,
#         k = 1.345,
#         tol = 1e-6,
#         maxIter = 10000
#     )
#
#     re <- saedevel:::optimizeRE.MSRFH(interfaceList)$fitre$x
#     out <- as.numeric(object$xy$x %*% object$beta + re)
#     attr(out, "re") <- as.numeric(re)
#     out
#
# }

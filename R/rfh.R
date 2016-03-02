#' Robust Fay Herriot Model
#'
#' @param formula (formula)
#' @param data (data.frame)
#' @param samplingVar (character)
#' @param ... arguments passed to methods
#' @param x0 (numeric) starting values
#' @param k (numeric) tuning constant
#' @param tol (numeric) numerical toloerance to be used during optimisation
#' @param y (numeric) response vector
#' @param X ([m|M]atrix) the design matrix
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
        public = c("call"),
        super = fitrfh(xy$y, xy$x, samplingVar, ...)
    )

}

#' @rdname rfh
#' @export
fitrfh <- function(y, X, samplingVar, x0 = c(rep(1, ncol(X)), 1), k = 1.345, tol = 1e-6) {
    # Non interactive fitting function for robust FH
    # y: (numeric) response
    # X: ((M|m)atrix) design matrix
    # samplingVar: (numeric) the sampling variances
    # x0: (numeric) starting values
    # k: (numeric) tuning constant

    # robust settings:
    psi <- . %>% psiOne(k)
    K <- getK(k)

    # Algorithm settings
    convCrit <- convCritRelative(tol)

    oneIter <- function(param) {
        param <- as.numeric(param)
        beta <- param[-length(param)]
        sigma2 <- param[length(param)]

        fpBeta <- fixedPointRobustBeta(
            y, X, matVFH(sigma2, samplingVar), psi = psi
        ) %>% addHistory

        beta <- fixedPoint(fpBeta, beta, convCrit)

        fpSigma2 <- fixedPointRobustDelta(
            y, X, beta, . %>% matVFH(samplingVar), psi, K
        ) %>% addConstraintMin(0) %>% addHistory

        sigma2 <- fixedPoint(fpSigma2, sigma2, convCrit)

        list(beta, sigma2)
    }

    # Fitting Model Parameter:
    out <- fixedPoint(addStorage(oneIter), x0, convCrit)

    beta <- out[-length(out)]
    variance <- out[length(out)]
    xy <- list(y = y, x = X)
    optimisationProgress <- storage$reformat(attr(out, "storage"))

    stripSelf(retList(
        "rfh",
        c("beta", "variance", "psi", "samplingVar", "xy", "optimisationProgress",
          "k", "tol", "K")))

}

#' @param object (rfh) an object of class rfh
#' @param type (character) one in \code{c("linear", "REBLUP")}
#' @param mse (character) which type of mse you want to compute for the
#'   predictions. See \link{mse} for available options.
#'
#' @rdname rfh
#' @export
predict.rfh <- function(object, type = "REBLUP", mse = "none", ...) {

    addPseudo <- function(out, isTrue, object, re, V) {
        # out: the data.frame with predictions
        # isTrue: (logical(1))
        if (isTrue) cbind(out, mse(object, "pseudo", re, V)["pseudo"])
        else out
    }

    addBoot <- function(out, isTrue, object, re, V, ...) {
        if (isTrue) cbind(out, mse(object, "boot", re, V, ...)["boot"])
        else out
    }

    if ("linear" == type) {
        mse <- "none" # no mse then
        re <- 0
    }

    if ("REBLUP" == type) {
        V <- variance(object)
        re <- fitReCCST(
            object$xy$y, object$xy$x, object$beta,
            variance(object)
        )
    }

    Xb <- object$xy$x %*% object$beta
    out <- data.frame(prediction = as.numeric(Xb + re))
    names(out) <- type
    out$re <- re
    out <- addPseudo(out, "pseudo" %in% mse, object, re, V)
    out <- addBoot(out, "boot" %in% mse, object, re, V, ...)

    out

}

fitReCCST <- function(y, X, beta, matV, psi = psiOne, convCrit = convCritAbsolute()) {
    # y: (numeric) response
    # X: (Matrix) Design-Matrix
    # beta: (numeric) estimated fixed effects
    # matV: (list) list of functions see e.g. matVFH
    # psi: (function) influence function

    # Non-robust random effects as starting values:
    startingValues <- as.numeric(tcrossprod(matV$Vu(), matV$Z()) %*% matV$VInv() %*%
                                     (y - X %*% beta))

    fpFun <- fixedPointRobustRandomEffect(
        y, X, beta, matV, psi
    )

    fixedPoint(fpFun, startingValues, convCrit)
}

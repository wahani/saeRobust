#' Model Variance Covariance Structures
#'
#' Functions to create the variance-covariance structure for the implemented
#' models. All these functions return a list with methods to compute variations.
#'
#' @param .sigma2 (numeric) a scalar. The variance parameter of the random
#'   effects.
#' @param .samplingVar (numeric) the vector of sampling variances of the direct
#'   estimator.
#' @param .rho (numeric) the correlation parameter for correlated random effects.
#' @param .W (matrix) proximity matrix.
#'
#' @details
#'
#' \code{matVFH} constructs variance-covariance matrix for a FH-model.
#'
#' \code{matVSFH} the spatial FH model.
#'
#' @rdname matModels
#' @export
matVFH <- function(.sigma2, .samplingVar) {

  .diag <- function(x) Diagonal(x = x)

  Vu <- getter(rep(.sigma2, length(.samplingVar)), .diag)
  VuInv <- getter(1 / rep(.sigma2, length(.samplingVar)), .diag)
  Ve <- getter(.samplingVar, .diag)
  VeInv <- getter(1 / .samplingVar, .diag)
  V <- getter(Vu() + Ve())
  VInv <- getter(solve(V()))
  Z <- getter(Diagonal(length(.samplingVar)))

  deriv <- list(
    getter(Diagonal(length(.samplingVar)))
  )

  retList()

}

#' @export
#' @rdname matModels
matVSFH <- function(.rho, .sigma2, .W, .samplingVar, .deriv = c("rho", "sigma2")) {

  .diag <- function(x) Diagonal(x = x)

  Ve <- getter(.samplingVar, .diag)
  VeInv <- getter(1 / .samplingVar, .diag)

  Omega1 <- getter(matOmega1(W = .W, rho = .rho), Matrix)
  Vu <- getter(.sigma2 * Omega1())
  VuInv <- getter(solve(Vu()))

  V <- getter(Vu() + Ve())
  VInv <- getter(solve(V()))

  Z <- getter(Diagonal(length(.samplingVar)))

  deriv <- list(
    rho = getter(matVDerR1(
      .rho, .sigma2, as.matrix(Z()), as.matrix(Omega1()), as.matrix(.W))),
    sigma2 = getter(Omega1())
  )[.deriv]

  retList()

}

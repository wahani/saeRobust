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
matVSFH <- function(.rho, .sigma2, .W, .samplingVar) {

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
  )

  retList()

}

#' @export
#' @rdname matModels
matVTFH <- function(.rho, .sigma2, .nTime, .samplingVar) {

  .diag <- function(x) Diagonal(x = x)
  .emptyMatrix <- function(dim) Matrix(0, dim, dim)
  .nDomains <- length(.samplingVar) / .nTime


  Ve <- getter(.samplingVar, .diag)
  VeInv <- getter(1 / .samplingVar, .diag)

  Omega1 <- getter(Diagonal(.nDomains))
  Omega2 <- getter(matOmega2(.nTime, .rho))

  Z <- getter(matTZ(.nDomains, .nTime))
  Z1 <- getter(matTZ1(.nDomains, .nTime))

  Vu <- getter(
    bdiag(.sigma2[1] * Omega1(),
          .sigma2[2] * matBlockDiagonal(Omega2(), .nDomains))
  )

  .V <- getter(matVInvT(
    as.matrix(Omega1()),
    .sigma2[1],
    .rho, .sigma2[2],
    as.matrix(Z1()),
    .samplingVar
  ))

  V <- getter(.V()$V)
  VInv <- getter(.V()$VInv)

  deriv <- list(
    sigma21 = getter(matVDerS1(as.matrix(Omega1()), as.matrix(Z1()))),
    rho = getter(matVDerR2(.rho, .sigma2[2], Omega2(), .nDomains)),
    sigma22 = getter(matVDerS2(Omega2(), .nDomains))
  )

  .ZVuZ <- getter(matVInvT(
    as.matrix(Omega1()),
    .sigma2[1],
    .rho, .sigma2[2],
    as.matrix(Z1()),
    rep(0, length(.samplingVar))
  ))

  ZVuZInv <- getter(.ZVuZ()$VInv)

  ZVuBarZ <- list(
    sigma21 = getter(Z() %*% tcrossprod(
      bdiag(Omega1(), .emptyMatrix(.nDomains * .nTime)), Z())),
    sigma22 = getter(Z() %*% tcrossprod(
      bdiag(.emptyMatrix(.nDomains), matBlockDiagonal(Omega2(), .nDomains)), Z()))
  )

  retList()

}

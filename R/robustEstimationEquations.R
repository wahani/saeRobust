robEstEqu <- function(y, X, .beta, u, matV, psi, K) {

  U <- getter(matU(matV$V()))
  Ue <- getter(matU(matV$Ve()))
  Uu <- getter(matU(matV$Vu()))

  psiResid <- getter(psi(U()$sqrtInv() %*% (y - X %*% .beta)))
  psiResidE <- getter(psi(Ue()$sqrtInv() %*% (y - X %*% .beta - matV$Z() %*% u)))
  psiResidU <- getter(psi(Uu()$sqrtInv() %*% u))

  beta <- getter({
    crossprod(X, matV$VInv()) %*% U()$sqrt() %*% psiResid()
  }, as.numeric)

  delta <- getter({
    lapply(
      matV$deriv,
      function(deriv) {
        browser()
        as.numeric(
          crossprod(psiResid(), U()$sqrt()) %*% matV$VInv() %*% deriv() %*%
            matV$VInv() %*% U()$sqrt() %*% psiResid() -
            matTrace(K * matV$VInv() %*% deriv())
        )
      })
  }, unlist)

  re <- getter({
    crossprod(matV$Z(), matV$VeInv()) %*% Ue()$sqrt() %*% psiResidE() -
      matV$VuInv() %*% Uu()$sqrt() %*% psiResidU()
  }, as.numeric)

  retList(public = c("beta", "delta", "re"))
}

#' Compute values of robust score functions
#'
#' These values should be close to zero for a model fit.
#'
#' @param object a fitted object
#' @param filter (character) a selection of values to be computed
#' @param ... arguments passed to methods
#'
#' @export
#' @rdname score
score <- function(object, filter, ...) UseMethod("score")

#' @export
score.default <- function(object, filter = c("beta", "delta", "re"), ...) {
  scores <- robEstEqu(
    object$y,
    object$x,
    object$coefficients,
    object$re,
    variance(object),
    object$psi,
    object$K)
  out <- lapply(filter, function(f) scores[[f]]())
  names(out) <- filter
  out
}

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

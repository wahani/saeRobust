robEstEqu <- function(y, X, .beta, u, matV, psi) {

    U <- getter(matU(matV$V()))
    Ue <- getter(matU(matV$Ve()))
    Uu <- getter(matU(matV$Vu()))

    psiResid <- getter(psi(U()$sqrtInv() %*% (y - X %*% .beta)))
    psiResidE <- getter(psi(Ue()$sqrtInv() %*% (y - X %*% .beta - u)))
    psiResidU <- getter(psi(Uu()$sqrtInv() %*% u))

    beta <- getter(
        crossprod(X, matV$VInv()) %*% U()$sqrt() %*% psiResid()
    )

    delta <- getter(
        crossprod(psiResid(), U()$sqrt()) %*% matV$VInv() %*% deriv %*%
            matV$VInv() %*% U()$sqrt() %*% psiResid() -
            matTrace(K %*% matV$VInv() %*% deriv)
    )

    re <- getter({
        crossprod(matV$Z(), matV$VeInv()) %*% Ue()$sqrt() %*% psiResidE() -
            matV$VuInv() %*% Uu()$sqrt() %*% psiResidU()
    })

    retList(public = c("beta", "delta", "re"))
}

set.seed(1)
X <- cbind(1, 1:10)
u <- rnorm(10)
y <- X %*% c(1, 1) + u
matV <- matVFH(var(u), rep(0.01, 10))

ree <- robEstEqu(y, X, c(1, 1), u, matV, identity)
ree$beta()
ree$delta()
ree$re()

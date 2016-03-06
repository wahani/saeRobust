context("Model Variance Structures")

test_that("FH", {

  expectEqual <- function(a, b) {
    testthat::expect_equal(a, b)
  }

  matV <- matVFH(0.5, c(1, 1))

  expectEqual(Matrix::Diagonal(x = 1 / (0.5 + c(1, 1))), matV$VInv())

})

test_that("SFH", {

  expectEqual <- function(a, b) {
    testthat::expect_equal(a, b, check.attributes = FALSE)
  }

  W <- spdep::nb2mat(spdep::cell2nb(3, 1, "rook"), style = "W")

  matV <- matVSFH(0.5, 0.5, W, c(1, 1, 1))
  Omega1 <- Matrix(
    solve(t(diag(3) - 0.5 * W) %*% (diag(3) - 0.5 * W)),
    forceCheck = TRUE
  )

  expectEqual(
    matV$VInv(),
    solve(0.5 * Omega1 + Diagonal(3))
  )

  expectEqual(
    as.matrix(matV$deriv$rho()),
    as.matrix(- 0.5 * Omega1 %*% (- W - t(W) + 2 * 0.5 * t(W) %*% W) %*% Omega1)
  )

  expectEqual(
    as.matrix(matV$deriv$sigma2()),
    as.matrix(Omega1)
  )

})

test_that("TFH", {

  expectEqual <- function(a, b) {
    testthat::expect_equal(a, b, check.attributes = FALSE)
  }

  matV <- matVTFH(0.5, c(2, 3), 3, rep(1, 9))
  Omega2 <- saeRobust:::matOmega2(3, 0.5)

  expectEqual(
    dim(matV$VInv()),
    c(9, 9)
  )

  expectEqual(
    as.matrix(matV$deriv$sigma21()),
    as.matrix(matTZ1(3, 3) %*% Diagonal(3) %*% t(matTZ1(3, 3)))
  )

  expectEqual(
    as.matrix(matV$deriv$rho()),
    saeRobust:::matVDerR2(0.5, 3, Omega2, 3)
  )

  expectEqual(
    as.matrix(matV$deriv$sigma22()),
    saeRobust:::matVDerS2(Omega2, 3)
  )

})



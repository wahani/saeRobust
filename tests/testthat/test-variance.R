context("variance")
test_that("matW", {

  expectEqual <- function(x, y) {
    testthat::expect_equal(x, y, tolerance = 1e-06)
  }

  library(saeSim)

  set.seed(1)
  dat <- base_id(10, 1) %>%
    sim_gen_e() %>%
    sim_gen_x() %>%
    sim_gen_v() %>%
    sim_resp_eq(y = 100 + 2 * x + v + e) %>%
    as.data.frame

  dat$samplingVar <- 1

  modelFit <- rfh(y ~ x, dat, "samplingVar")
  out <- predict(modelFit)

  # This should be the same as Xb + u which is used in predict:
  W <-  weights(modelFit)$W

  prediction <- W %*% modelFit$y

  expectEqual(
    as.numeric(out$REBLUP),
    as.numeric(prediction)
  )

  expectEqual(
    as.numeric(rowSums(W)),
    rep(1, 10)
  )

  expectEqual(
    modelFit$x %*% weights(modelFit)$A + # XA
      weights(modelFit)$B %*% ( # B
        Diagonal(length(modelFit$y)) - # I -
          modelFit$x %*% weights(modelFit)$A), # XA
    weights(modelFit)$W
  )

})

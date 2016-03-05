context("fit")

test_that("rfhfit is working", {
  library("saeSim")
  set.seed(1)
  dat <- base_id(10, 1) %>%
    sim_gen_e() %>%
    sim_gen_x() %>%
    sim_gen_v() %>%
    sim_resp_eq(y = 100 + 2 * x + v + e) %>%
    as.data.frame

  y <- dat$y
  X <- cbind(1, dat$x)
  samplingVar <- rep(16, nrow(dat))

  out <- saeRobust:::fitrfh(y, X, samplingVar)
  expect_is(out, "list")
  expect_is(out$coefficients, "numeric")
  expect_is(out$variance, "numeric")
})

test_that("fitrsfh is working", {
  library("saeSim")
  set.seed(2)
  nDomains <- 40
  dat <- base_id(nDomains, 1) %>%
    sim_gen_e() %>%
    sim_gen_x() %>%
    sim_gen(gen_v_sar(rho = 0.5, sd = 4, name = "v")) %>%
    sim_resp_eq(y = 100 + 2 * x + v + e) %>%
    as.data.frame

  y <- dat$y
  X <- cbind(1, dat$x)
  samplingVar <- rep(16, nrow(dat))
  W <- spdep::nb2mat(spdep::cell2nb(nDomains, 1, "rook"), style = "W")

  out <- saeRobust:::fitrsfh(y, X, samplingVar, W, x0Var = c(0.01, 1))
  # out$variance
  # out$iterations$correlation
  # score(out)$delta

  expect_is(out, "list")
  expect_is(out$coefficients, "numeric")
  expect_is(out$variance, "numeric")

})

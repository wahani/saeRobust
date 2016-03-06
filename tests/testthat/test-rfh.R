context("rfh")

# Note that I do not test the correctness of the estimators. Here I only check
# if the public representation of the model object looks how I want them to. You
# can find some tests on the correctness of the estimation equations.

test_that("rfh is working", {

    library("saeSim")
    set.seed(1)
    dat <- base_id(10, 1) %>%
        sim_gen_e() %>%
        sim_gen_x() %>%
        sim_gen_v() %>%
        sim_resp_eq(y = 100 + 2 * x + v + e) %>%
        as.data.frame

    dat$samplingVar <- 1

    out <- rfh(y ~ x, dat, "samplingVar", maxIter = 1, maxIterRe = 1)
    expect_is(out, "list")
    expect_is(out$coefficients, "numeric")
    expect_is(out$variance, "numeric")
    expect_is(out$samplingVar, "numeric")
    expect_is(out$y, "numeric")

})

test_that("spatial rfh is working", {

  library("saeSim")
  set.seed(1)
  dat <- base_id(10, 1) %>%
    sim_gen_e() %>%
    sim_gen_x() %>%
    sim_gen_v() %>%
    sim_resp_eq(y = 100 + 2 * x + v + e) %>%
    as.data.frame

  dat$samplingVar <- 1

  W <- testRook(10)

  out <- rfh(y ~ x, dat, "samplingVar", corSAR1(W), maxIter = 1, maxIterRe = 1)
  expect_is(out, "list")
  expect_is(out$coefficients, "numeric")
  expect_is(out$variance, "numeric")
  expect_is(out$samplingVar, "numeric")
  expect_is(out$y, "numeric")
  testthat::expect_is(out$W, "matrix")

})


test_that("predict.rfh", {

    expectIs <- function(x, a) {
        testthat::expect_is(x, a)
    }

    expectEqual <- function(x, y) {
        testthat::expect_equal(x, y)
    }

    library("saeSim")
    set.seed(1)
    dat <- base_id(10, 1) %>%
        sim_gen_e(sd = 1) %>%
        sim_gen_x() %>%
        sim_gen_v(sd = 1) %>%
        sim_resp_eq(y = 100 + 2 * x + v + e) %>%
        as.data.frame

    dat$samplingVar <- 1

    modelFit <- rfh(y ~ x, dat, "samplingVar")
    out <- predict(modelFit)

    expectIs(out, "data.frame")
    expectIs(out$re, "numeric")
    expectEqual(
        as.numeric(out$REBLUP - out$re),
        as.numeric(modelFit$x %*% modelFit$coefficients)
    )

    out <- predict(modelFit)
    expectEqual(names(out), c("REBLUP", "re"))

    out <- predict(modelFit, "linear")
    expectEqual(names(out), c("linear", "re"))

    out <- predict(modelFit)
    expectEqual(names(out), c("REBLUP", "re"))

})



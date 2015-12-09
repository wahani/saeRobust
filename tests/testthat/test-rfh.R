context("rfh")

# Note that I do not test the correctness of the estimators. I rely on the tests
# implementing for the objective functions and the algorithms.

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

    out <- saeRobustTools:::fitrfh(y, X, samplingVar)
    expect_is(out, "list")
    expect_is(out$beta, "numeric")
    expect_is(out$variance, "numeric")
})


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

    out <- rfh(y ~ x, dat, "samplingVar")
    expect_is(out, "list")
    expect_is(out$beta, "numeric")
    expect_is(out$variance, "numeric")
    expect_is(out$samplingVar, "numeric")
    expect_is(out$xy, "list")

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
        sim_gen_e() %>%
        sim_gen_x() %>%
        sim_gen_v() %>%
        sim_resp_eq(y = 100 + 2 * x + v + e) %>%
        as.data.frame

    dat$samplingVar <- 1

    modelFit <- rfh(y ~ x, dat, "samplingVar")
    out <- predict(modelFit)

    expectIs(out, "numeric")
    expectIs(attr(out, "re"), "numeric")
    expectEqual(
        as.numeric(out - attr(out, "re")),
        as.numeric(modelFit$xy$x %*% modelFit$beta)
    )
})

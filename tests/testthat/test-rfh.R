context("rfh")
test_that("rfh is working", {
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

    out <- rfh(y, X, samplingVar)
    expect_is(out, "list")
    expect_is(out$beta, "numeric")
    expect_is(out$variance, "numeric")
})

context("mat")
test_that("matW", {

    expectEqual <- function(x, y) {
        testthat::expect_equal(x, y)
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
    prediction <-
        matW(
            y = modelFit$xy$y,
            X = modelFit$xy$x,
            beta = modelFit$beta,
            u = attr(out, "re"),
            Diagonal(10, x = modelFit$variance + 1),
            VuSqrtInv = sqrt(solve(Diagonal(10, modelFit$variance))),
            VeSqrtInv = Diagonal(10, 1),
            psi = psiOne
        ) %*% modelFit$xy$y

    expectEqual(
        as.numeric(out),
        as.numeric(prediction)
    )

})

# devtools::install_github("wahani/dat")
library(dat)
library(saeSim)
library(saeRobustTools)
library(Matrix)
library(ggplot2)


comp_rfh <- function(dat) {
    modelFit <- rfh(y ~ x, dat, "samplingVar")
    prediction <- predict(modelFit)
    dat$rfh <- as.numeric(prediction)

    W <-  matW(
        y = modelFit$xy$y,
        X = modelFit$xy$x,
        beta = modelFit$beta,
        u = attr(prediction, "re"),
        Diagonal(length(modelFit$xy$y), x = modelFit$variance + 1),
        VuSqrtInv = sqrt(solve(Diagonal(length(modelFit$xy$y), modelFit$variance))),
        VeSqrtInv = sqrt(solve(Diagonal(length(modelFit$xy$y), modelFit$samplingVar))),
        psi = psiOne
    )

    dat$mse <- as.numeric(W^2 %*% modelFit$samplingVar + (dat$rfh - dat$y)^2)
    dat

}

setup <- base_id(40, 1) %>%
    sim_gen_x() %>%
    sim_gen_v(sd = 4) %>%
    as.data.frame %>%
    sim_gen_e() %>%
    sim_comp_pop(comp_var(samplingVar = 16)) %>%
    sim_resp_eq(y = 100 + 2 * x + v + e, trueVal = 100 + 2 * x + v) %>%
    sim_comp_agg(comp_rfh)

simDat <- setup %>% sim(100) %>% reduce(Dots(rbind))


ggDat <- simDat %>%
    mutar(MCMSE ~ mean((y - trueVal)^2),
          MSE ~ median(mse),
          by = "idD")

ggDat <- reshape2::melt(ggDat, id.vars = c("idD"))

ggplot(ggDat, aes(x = idD, y = value, colour = variable)) +
    geom_line()

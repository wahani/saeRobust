#' Plots
#'
#' Various implementations of diagnostic plots. They are linked together using the \link{plot} generic function.
#'
#' @param sample (numeric) a vector
#' @param x an object
#' @param y for mse estimates a filter for the predictors; otherwise ignored
#' @param alpha (numeric) between 0 and 1 - used in computation of confidence
#'   interval
#' @param ... ignored
#'
#' @details
#'
#' \code{qqPlot} a QQ-Plot using ggplot2
#'
#' \code{blandAltmanPlot} a Bland-Altman plot. Solid line is the mean. Dashed
#' lines are the upper and lower bound of a confidence interval assuming a
#' normal distribution. The alpha level can be set using \code{alpha}.
#'
#' @export
#' @rdname plot
#'
#' @examples
#' qqPlot(rnorm(10))
#' blandAltmanPlot(rnorm(10), rnorm(10))
plot.rfh <- function(x, y, ...) {
  # QQ: Sampling Error
  ggSamplingError <-
    qqPlot(x$residuals / sqrt(x$samplingVar)) +
    labs(
      y = "residuals / sqrt(samplingVar)",
      title = "Quantile-Quantile Plot for Standardised Residuals"
    )

  # QQ: Random Effects
  ggRandomEffects <-
    qqPlot(x$reblup - x$fitted) + # because 're' can have different length
    labs(
      y = "Random Effects",
      title = "Quantile-Quantile Plot for Random Effects"
    )

  #
  ggPlots <- list(
    residuals = ggSamplingError,
    randomEffects = ggRandomEffects
  )

  if (interactive()) {
    lapply(ggPlots, function(gg) {
      print(gg)
      dump <- readline("Press Enter for next plot.")
    })
  }

  invisible(ggPlots)

}

#' @method plot prediction.fitrfh
#' @rdname plot
#' @export
plot.prediction.fitrfh <- function(x, y, alpha = 0.05, ...) {

  varsToPlot <- names(x)[names(x) %in% c("linear", "reblup", "reblupbc")]
  ggPlots <- lapply(varsToPlot, function(varName) {
    xLab <- paste0("(direct + ", varName, ") / 2")
    yLab <- paste0("direct - ", varName)
    main <- paste0("Bland-Altman Plot: direct vs. ", varName)
    blandAltmanPlot(x[["direct"]], x[[varName]]) +
      labs(x = xLab, y = yLab, title = main)
  })

  if (interactive()) {
    lapply(ggPlots, function(gg) {
      print(gg)
      dump <- readline("Press Enter for next plot.")
    })
  }

  invisible(ggPlots)

}

#' @method plot mse.fitrfh
#' @rdname plot
#' @export
plot.mse.fitrfh <- function(x, y = "pseudo", ...) {
  # cv: coefficient of variation -- se / mean

  varToPlot <- names(x)[names(x) %in% c("reblup", "reblupbc")]
  mseToPlot <- names(x)[names(x) %in% paste0(y, c("", "bc"))]
  x <- x[order(x$samplingVar), ]

  ggDat <- data.frame(
    predictions = c(x$direct, unlist(as.list(x[varToPlot]))),
    mse = c(x$samplingVar, unlist(as.list(x[mseToPlot]))),
    method = c(rep("direct", NROW(x)), rep(varToPlot, each = NROW(x))),
    area = 1:NROW(x)
  )

  ggDat$cv <- 100 * abs(sqrt(ggDat$mse) / ggDat$predictions)

  cvPlot <- ggplot(ggDat, aes_string(x = "area", y = "cv", colour = "method")) +
    geom_point() +
    labs(
      x = "domain (sorted by increasing sampling variance of direct estimator)",
      y = "CV",
      title = c("Scatterplot for CV against sorted domain")
    )

  boxPlot <- ggplot(ggDat, aes_string(x = "method", y = "cv")) +
    geom_boxplot() +
    coord_flip() +
    labs(x = NULL, y = "CV", title = "Boxplot for CV")

  ggPlots <- list(scatterPlot = cvPlot, boxPlot = boxPlot)

  if (interactive()) {
    lapply(ggPlots, function(gg) {
      print(gg)
      dump <- readline("Press Enter for next plot.")
    })
  }

  invisible(ggPlots)

}


#' @export
#' @rdname plot
qqPlot <- function(sample) {

  y <- quantile(sample[!is.na(sample)], c(0.25, 0.75))
  x <- qnorm(c(0.25, 0.75))
  slope <- diff(y)/diff(x)
  int <- y[1L] - slope * x[1L]

  ggDat <- data.frame(sample = sample)

  ggplot(ggDat, aes(sample = sample)) +
    geom_qq() +
    geom_abline(slope = slope, intercept = int) +
    labs(title = "Quantile-Quantile Plot")

}

#' @export
#' @rdname plot
blandAltmanPlot <- function(x, y, alpha = 0.05) {

  ggDat <- data.frame(x = x, y = y)
  me <- mean((x - y), na.rm = TRUE)
  width <- qnorm(alpha / 2) * sd(x - y, TRUE)
  lower <- me - width
  upper <- me + width

  ggplot(ggDat, aes(x = (x + y) / 2, y = x - y))+
    geom_point() +
    geom_hline(yintercept = me, linetype = 1) +
    geom_hline(yintercept = lower, linetype = 2) +
    geom_hline(yintercept = upper, linetype = 2) +
    labs(title = "Bland-Altman Plot")

}

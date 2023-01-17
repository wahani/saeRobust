library("saeRobust")

if (requireNamespace("testthat", quietly = TRUE))
  if (requireNamespace("sae", quietly = TRUE))
    if (requireNamespace("saeSim", quietly = TRUE))
      testthat::test_check("saeRobust")

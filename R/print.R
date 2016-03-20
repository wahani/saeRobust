#' @method print rfh
#' @export
print.rfh <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  # Call:
  cat(
    "\nCall:\n",
    paste(deparse(x$call), sep = "\n", collapse = "\n"),
    "\n\n", sep = ""
  )

  # Coefficients:
  if (length(coef(x))) {
    cat("Coefficients:\n")
    print.default(
      format(coef(x), digits = digits),
      print.gap = 2L, quote = FALSE
    )
  }
  else cat("No coefficients\n")
  cat("\n")

  # Variance Components:
  cat("Variance Components:\n")
  print.default(
    format(x$variance, digits = digits),
    print.gap = 2L, quote = FALSE
  )
  cat("\n")

  invisible(x)

}

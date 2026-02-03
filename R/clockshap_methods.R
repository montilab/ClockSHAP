#' @export
print.clockshap <- function(x, ...) {

  cat("<ClockSHAP object>\n")
  cat("  Samples :", nrow(x$phi), "\n")
  cat("  Features:", ncol(x$phi), "\n")
  cat("\n")

  dev_range <- range(x$deviation)
  cat("  Deviation range:\n")
  cat("    Min:", round(dev_range[1], 3), "\n")
  cat("    Max:", round(dev_range[2], 3), "\n")

  invisible(x)
}

#' @export
summary.clockshap <- function(object, ...) {

  list(
    n_samples  = nrow(object$phi),
    n_features = ncol(object$phi),
    deviation  = summary(object$deviation),
    age        = summary(object$age)
  )
}

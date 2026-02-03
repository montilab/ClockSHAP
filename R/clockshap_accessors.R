#' Extract ClockSHAP feature contributions
#'
#' @param x A `clockshap` object.
#'
#' @return A numeric matrix of feature contributions
#'   (samples x features).
#' @export
clockshap_phi <- function(x) {
  if (!inherits(x, "clockshap")) {
    stop("`x` must be a clockshap object.", call. = FALSE)
  }
  x$phi
}


#' Extract ClockSHAP deviation values
#'
#' @param x A `clockshap` object.
#'
#' @return A numeric vector of deviations.
#' @export
clockshap_deviation <- function(x) {
  if (!inherits(x, "clockshap")) {
    stop("`x` must be a clockshap object.", call. = FALSE)
  }
  x$deviation
}


#' Extract ClockSHAP predicted ages
#'
#' @param x A `clockshap` object.
#'
#' @return A numeric vector of predicted ages.
#' @export
clockshap_predicted <- function(x) {
  if (!inherits(x, "clockshap")) {
    stop("`x` must be a clockshap object.", call. = FALSE)
  }
  x$predicted
}


#' Extract ClockSHAP expected ages
#'
#' @param x A `clockshap` object.
#'
#' @return A numeric vector of expected ages.
#' @export
clockshap_expected <- function(x) {
  if (!inherits(x, "clockshap")) {
    stop("`x` must be a clockshap object.", call. = FALSE)
  }
  x$expected
}

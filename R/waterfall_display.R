#' Compute reconciled values for a ClockSHAP waterfall
#'
#' Prepare rounded values for waterfall annotations while preserving the two
#' displayed identities: feature contributions sum to the deviation, and the
#' expected age plus the deviation equals the predicted age.
#'
#' @param x A `clockshap` object.
#' @param sample A sample index (integer) or row name (character).
#' @param top_n Number of top absolute contributions to display.
#' @param digits Number of decimal places used for all numeric labels.
#'
#' @return A list containing the displayed `expected`, `predicted`, `deviation`,
#'   and named `contributions` values.
#'
#' @details To preserve exact arithmetic after rounding, an endpoint or
#'   contribution label may differ from independently rounding its stored value
#'   by one unit in the last displayed decimal place.
#'
#' @export
clockshap_waterfall_values <- function(x, sample, top_n = 10, digits = 1) {
  if (!inherits(x, "clockshap")) {
    stop("`x` must be a clockshap object.", call. = FALSE)
  }
  top_n <- .validate_waterfall_top_n(top_n)
  digits <- .validate_waterfall_digits(digits)

  phi <- clockshap_phi(x)
  i <- .resolve_waterfall_sample(phi, sample)
  contributions <- .waterfall_contributions(phi[i, ], top_n)

  expected <- unname(clockshap_expected(x)[i])
  predicted <- unname(clockshap_predicted(x)[i])
  deviation <- unname(clockshap_deviation(x)[i])
  if (!all(is.finite(c(expected, predicted, deviation, contributions)))) {
    stop("Selected sample has non-finite waterfall values; cannot format.",
         call. = FALSE)
  }
  if (!isTRUE(all.equal(
    unname(expected + deviation), unname(predicted), tolerance = 1e-10
  ))) {
    stop("ClockSHAP invariant violated: expected plus deviation != predicted.",
         call. = FALSE)
  }
  if (!isTRUE(all.equal(
    unname(sum(contributions)), unname(deviation), tolerance = 1e-10
  ))) {
    stop("ClockSHAP invariant violated: contributions do not sum to deviation.",
         call. = FALSE)
  }

  scale <- 10^digits
  deviation_units <- round(deviation * scale)
  expected_units <- round(expected * scale)
  contribution_units <- .reconcile_contribution_units(
    contributions,
    deviation_units,
    scale
  )

  displayed <- list(
    expected = expected_units / scale,
    predicted = (expected_units + deviation_units) / scale,
    deviation = deviation_units / scale,
    contributions = contribution_units / scale
  )
  displayed$contributions[displayed$contributions == 0] <- 0
  displayed
}


.validate_waterfall_top_n <- function(top_n) {
  if (!is.numeric(top_n) || length(top_n) != 1 || !is.finite(top_n) ||
      top_n < 1 || top_n != as.integer(top_n)) {
    stop("`top_n` must be a positive integer.", call. = FALSE)
  }
  as.integer(top_n)
}


.validate_waterfall_digits <- function(digits, argument = "digits") {
  if (!is.numeric(digits) || length(digits) != 1 || !is.finite(digits) ||
      digits < 0 || digits > 9 || digits != as.integer(digits)) {
    stop(sprintf("`%s` must be an integer between 0 and 9.", argument),
         call. = FALSE)
  }
  as.integer(digits)
}


.resolve_waterfall_sample <- function(phi, sample) {
  if (is.character(sample)) {
    if (length(sample) != 1 || is.na(sample) || is.null(rownames(phi)) ||
        !sample %in% rownames(phi)) {
      stop("Sample name not found in ClockSHAP object.", call. = FALSE)
    }
    return(match(sample, rownames(phi)))
  }

  if (!is.numeric(sample) || length(sample) != 1 || !is.finite(sample) ||
      sample != as.integer(sample)) {
    stop("`sample` must be a single sample name or integer index.",
         call. = FALSE)
  }
  i <- as.integer(sample)
  if (i < 1 || i > nrow(phi)) {
    stop("Sample index out of bounds.", call. = FALSE)
  }
  i
}


.waterfall_contributions <- function(phi, top_n) {
  phi <- phi[order(-abs(phi))]
  top <- head(phi, top_n)
  other <- sum(phi) - sum(top)
  c(top, Other = other)
}


.reconcile_contribution_units <- function(values, target_units, scale) {
  raw_units <- unname(values) * scale
  units <- round(raw_units)
  residual <- as.integer(round(target_units - sum(units)))

  if (residual != 0L) {
    direction <- sign(residual)
    remainders <- raw_units - units
    feature_names <- names(values)
    if (is.null(feature_names)) {
      feature_names <- rep("", length(values))
    }
    order_index <- if (direction > 0) {
      order(-remainders, -abs(values), feature_names, seq_along(values))
    } else {
      order(remainders, -abs(values), feature_names, seq_along(values))
    }

    eligible <- vapply(order_index, function(j) {
      proposed <- units[j] + direction
      (values[j] > 0 && proposed >= 0) ||
        (values[j] < 0 && proposed <= 0) ||
        values[j] == 0
    }, logical(1))
    recipients <- order_index[eligible]

    if (length(recipients) < abs(residual)) {
      stop("Unable to reconcile contribution labels without changing a sign.",
           call. = FALSE)
    }
    recipients <- head(recipients, abs(residual))
    units[recipients] <- units[recipients] + direction
  }

  names(units) <- names(values)
  units[units == 0] <- 0
  units
}


.format_waterfall_signed <- function(values, digits) {
  labels <- sprintf(paste0("%+.", digits, "f"), values)
  labels[values == 0] <- sprintf(paste0("%.", digits, "f"), 0)
  labels
}

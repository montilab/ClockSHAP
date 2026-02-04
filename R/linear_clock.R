## Linear aging clock representation

.make_linear_clock <- function(alpha, beta, mu, sigma) {
  .validate_named_numeric(beta, "beta")
  .validate_named_numeric(mu, "mu")
  .validate_named_numeric(sigma, "sigma")

  ## require same name set (order can differ)
  .validate_same_names(beta, mu, sigma)

  if (!is.numeric(alpha) || length(alpha) != 1 || !is.finite(alpha)) {
    stop("`alpha` must be a single finite numeric intercept.", call. = FALSE)
  }

  ## align preprocessing vectors to beta order
  feat_names <- names(beta)
  mu    <- .reorder_named(mu, feat_names, name = "mu")
  sigma <- .reorder_named(sigma, feat_names, name = "sigma")

  if (any(sigma <= 0)) {
    stop("All elements of `sigma` must be positive (non-zero).", call. = FALSE)
  }

  clock <- list(
    alpha = alpha,
    beta  = beta,
    mu    = mu,
    sigma = sigma
  )

  class(clock) <- "linear_clock"
  clock
}

#' Create a linear aging clock
#'
#' Construct a linear aging clock specification from coefficients and
#' preprocessing parameters. This does not train a model.
#'
#' @param alpha Numeric scalar intercept.
#' @param beta Named numeric vector of feature coefficients.
#' @param mu Named numeric vector of training-set feature means.
#' @param sigma Named numeric vector of training-set feature standard deviations.
#'
#' @return An object of class `linear_clock`.
#' @export
linear_clock <- function(alpha, beta, mu, sigma) {
  .make_linear_clock(alpha = alpha, beta = beta, mu = mu, sigma = sigma)
}

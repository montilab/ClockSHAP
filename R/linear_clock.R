## Linear aging clock representation

.make_linear_clock <- function(alpha, beta, mu, sigma) {
  .validate_named_numeric(beta, "beta")
  .validate_named_numeric(mu, "mu")
  .validate_named_numeric(sigma, "sigma")

  .validate_same_names(beta, mu, sigma)

  if (!is.numeric(alpha) || length(alpha) != 1) {
    stop("`alpha` must be a single numeric intercept.", call. = FALSE)
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

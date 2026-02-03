#' ClockSHAP decomposition
#'
#' Decompose biological age deviation into exact, age-conditioned feature
#' contributions for a linear aging clock.
#'
#' @param features A matrix or data.frame of feature values
#'   (samples x features).
#' @param age Numeric vector of chronological ages (length = nrow(features)).
#' @param clock A `linear_clock` object.
#' @param reference A `reference_profile` object.
#'
#' @return An object of class `clockshap` containing per-feature contributions,
#'   deviation values, and metadata.
#'
#' @examples
#' set.seed(1)
#'
#' features <- matrix(rnorm(20), nrow = 5, ncol = 4)
#' colnames(features) <- paste0("F", 1:4)
#' age <- c(40, 50, 60, 70, 80)
#'
#' clock <- linear_clock(
#'   alpha = 10,
#'   beta  = setNames(runif(4), colnames(features)),
#'   mu    = setNames(rep(0, 4), colnames(features)),
#'   sigma = setNames(rep(1, 4), colnames(features))
#' )
#'
#' ref <- fit_reference_profile(
#' features = as.data.frame(features),
#' age = age
#' )
#'
#' cs <- clockshap(features, age, clock, ref)
#' cs
#' summary(cs)
#'
#'
#' @export

clockshap <- function(features, age, clock, reference) {

  ## ---------------------- validation ----------------------

  if (!inherits(clock, "linear_clock")) {
    stop("`clock` must be a `linear_clock` object.", call. = FALSE)
  }
  if (!inherits(reference, "reference_profile")) {
    stop("`reference` must be a `reference_profile` object.", call. = FALSE)
  }

  if (!is.matrix(features) && !is.data.frame(features)) {
    stop("`features` must be a matrix or data.frame.", call. = FALSE)
  }

  features <- as.matrix(features)

  if (!is.numeric(age)) {
    stop("`age` must be numeric.", call. = FALSE)
  }
  if (nrow(features) != length(age)) {
    stop("Length of `age` must match number of rows in `features`.",
         call. = FALSE)
  }

  .validate_same_names(
    clock$beta,
    clock$mu,
    clock$sigma,
    reference$gamma0,
    reference$gamma1
  )

  if (!identical(names(clock$beta), colnames(features))) {
    stop("Feature names in `features` must match clock and reference names.",
         call. = FALSE)
  }


  ## ---------------------- unpack objects ----------------------

  beta  <- clock$beta
  mu    <- clock$mu
  sigma <- clock$sigma
  alpha <- clock$alpha

  gamma0 <- reference$gamma0
  gamma1 <- reference$gamma1

  ## ---------------------- standardize features ----------------------

  X_std <- sweep(features, 2, mu, "-")
  X_std <- sweep(X_std, 2, sigma, "/")

  ## ---------------------- expected features at age ----------------------

  X_exp <- sweep(
    matrix(gamma1, nrow = nrow(features), ncol = length(gamma1), byrow = TRUE),
    1,
    age,
    "*"
  )
  X_exp <- sweep(X_exp, 2, gamma0, "+")
  colnames(X_exp) <- names(gamma0)

  X_exp_std <- sweep(X_exp, 2, mu, "-")
  X_exp_std <- sweep(X_exp_std, 2, sigma, "/")

  ## ---------------------- ClockSHAP contributions ----------------------

  phi <- sweep(X_std - X_exp_std, 2, beta, "*")

  ## ---------------------- predictions and deviation ----------------------

  y_hat <- alpha + as.numeric(X_std %*% beta)
  y_exp <- alpha + as.numeric(X_exp_std %*% beta)
  delta <- y_hat - y_exp

  ## ---------------------- invariant check ----------------------

  sum_phi <- rowSums(phi)
  if (!isTRUE(all.equal(sum_phi, delta, tol = 1e-10))) {
    stop("ClockSHAP invariant violated: sum(phi) != deviation.",
         call. = FALSE)
  }

  ## ---------------------- return object ----------------------

  out <- list(
    phi       = phi,
    deviation = delta,
    predicted = y_hat,
    expected  = y_exp,
    age       = age
  )

  class(out) <- "clockshap"
  out
}

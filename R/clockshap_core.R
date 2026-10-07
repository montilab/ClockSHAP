#' ClockSHAP decomposition
#'
#' Decompose age deviation into exact, age-conditioned feature
#' contributions for a linear aging clock.
#'
#' @param features A matrix or data.frame of feature values
#'   (samples x features).
#' @param age Numeric vector of chronological ages (length = nrow(features)).
#'   If sample identifiers are provided (`rownames(features)` and
#'   `names(age)`), they must match exactly and in the same order.
#' @param clock A `linear_clock` object.
#' @param reference A `reference_profile` object.
#' @param max_age Optional numeric scalar. Samples with `age > max_age` are
#'   treated as out of reference range and return `NA` for `expected`,
#'   `deviation`, and `phi`.
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
#'   features = as.data.frame(features),
#'   age = age
#' )
#'
#' cs <- clockshap(features, age, clock, ref)
#' cs
#' summary(cs)
#'
#' @export
clockshap <- function(features, age, clock, reference, max_age = NULL) {

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

  if (is.null(colnames(features))) {
    stop("`features` must have column names matching the clock features.",
         call. = FALSE)
  }

  if (anyNA(features)) {
    stop("`features` must not contain missing values.", call. = FALSE)
  }

  if (!is.numeric(age)) {
    stop("`age` must be numeric.", call. = FALSE)
  }
  if (length(age) != nrow(features)) {
    stop("Length of `age` must match number of rows in `features`.",
         call. = FALSE)
  }
  if (anyNA(age)) {
    stop("`age` must not contain missing values.", call. = FALSE)
  }

  feat_rows <- rownames(features)
  age_names <- names(age)
  if (!is.null(feat_rows) || !is.null(age_names)) {
    if (is.null(feat_rows)) {
      stop("`rownames(features)` are required when `age` is named.",
           call. = FALSE)
    }
    if (is.null(age_names)) {
      stop("`names(age)` are required when `features` has row names.",
           call. = FALSE)
    }
    if (anyNA(feat_rows) || any(feat_rows == "") || anyDuplicated(feat_rows)) {
      stop("`rownames(features)` must be non-missing, non-empty, and unique.",
           call. = FALSE)
    }
    if (anyNA(age_names) || any(age_names == "") || anyDuplicated(age_names)) {
      stop("`names(age)` must be non-missing, non-empty, and unique.",
           call. = FALSE)
    }
    if (!identical(age_names, feat_rows)) {
      stop("`names(age)` must exactly match `rownames(features)` in the same order.",
           call. = FALSE)
    }
  }
  if (!is.null(max_age)) {
    if (!is.numeric(max_age) || length(max_age) != 1 || !is.finite(max_age)) {
      stop("`max_age` must be NULL or a single finite numeric value.",
           call. = FALSE)
    }
  }

  ## require same feature-name set across clock and reference
  .validate_same_names(
    clock$beta,
    clock$mu,
    clock$sigma,
    reference$gamma0,
    reference$gamma1
  )

  ## ---------------------- align feature order ----------------------

  feat_names <- names(clock$beta)

  ## require exact match between features and clock names (order can differ)
  if (!identical(sort(colnames(features)), sort(feat_names))) {
    missing <- setdiff(feat_names, colnames(features))
    extra   <- setdiff(colnames(features), feat_names)

    msg <- "Feature names in `features` must match clock and reference names."
    if (length(missing) > 0) {
      msg <- paste0(msg, " Missing: ", paste(missing, collapse = ", "), ".")
    }
    if (length(extra) > 0) {
      msg <- paste0(msg, " Extra: ", paste(extra, collapse = ", "), ".")
    }
    stop(msg, call. = FALSE)
  }

  ## reorder to the clock's canonical feature order
  features <- features[, feat_names, drop = FALSE]

  beta  <- .reorder_named(clock$beta,  feat_names, name = "clock$beta")
  mu    <- .reorder_named(clock$mu,    feat_names, name = "clock$mu")
  sigma <- .reorder_named(clock$sigma, feat_names, name = "clock$sigma")
  alpha <- clock$alpha

  gamma0 <- .reorder_named(reference$gamma0, feat_names, name = "reference$gamma0")
  gamma1 <- .reorder_named(reference$gamma1, feat_names, name = "reference$gamma1")

  ## ---------------------- standardize features ----------------------

  X_std <- sweep(features, 2, mu, "-")
  X_std <- sweep(X_std, 2, sigma, "/")

  ## ---------------------- expected features at age ----------------------
  valid_ref <- rep(TRUE, length(age))
  age_ref <- age
  if (!is.null(max_age)) {
    valid_ref <- age <= max_age
    age_ref <- pmin(age, max_age)
  }

  X_exp <- sweep(
    matrix(gamma1, nrow = nrow(features), ncol = length(gamma1), byrow = TRUE),
    1,
    age_ref,
    "*"
  )
  X_exp <- sweep(X_exp, 2, gamma0, "+")
  colnames(X_exp) <- feat_names

  X_exp_std <- sweep(X_exp, 2, mu, "-")
  X_exp_std <- sweep(X_exp_std, 2, sigma, "/")

  ## ---------------------- ClockSHAP contributions ----------------------

  phi <- sweep(X_std - X_exp_std, 2, beta, "*")

  ## ---------------------- predictions and deviation ----------------------

  y_hat <- alpha + as.numeric(X_std %*% beta)
  y_exp <- alpha + as.numeric(X_exp_std %*% beta)
  delta <- y_hat - y_exp

  if (!all(valid_ref)) {
    y_exp[!valid_ref] <- NA_real_
    delta[!valid_ref] <- NA_real_
    phi[!valid_ref, ] <- NA_real_
  }


  ## ---------------------- invariant check ----------------------

  if (any(valid_ref)) {
    sum_phi <- rowSums(phi[valid_ref, , drop = FALSE])
    if (!isTRUE(all.equal(
      unname(sum_phi),
      unname(delta[valid_ref]),
      tolerance = 1e-10
    ))) {
      stop("ClockSHAP invariant violated: sum(phi) != deviation.",
           call. = FALSE)
    }
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

## Reference age-feature relationships

.fit_reference_trends <- function(features, age) {
  if (!is.matrix(features) && !is.data.frame(features)) {
    stop("`features` must be a matrix or data.frame (samples x features).",
         call. = FALSE)
  }

  if (!is.numeric(age)) {
    stop("`age` must be numeric.", call. = FALSE)
  }

  features <- as.matrix(features)

  if (nrow(features) != length(age)) {
    stop("Number of samples in `features` must match length of `age`.",
         call. = FALSE)
  }

  if (anyNA(age)) {
    stop("`age` must not contain missing values.", call. = FALSE)
  }
  if (anyNA(features)) {
    stop("`features` must not contain missing values.", call. = FALSE)
  }

  ## Ensure feature names exist (required for alignment downstream)
  if (is.null(colnames(features))) {
    colnames(features) <- paste0("V", seq_len(ncol(features)))
  }

  gamma0 <- gamma1 <- numeric(ncol(features))

  for (k in seq_len(ncol(features))) {
    fit <- lm(features[, k] ~ age)
    gamma0[k] <- coef(fit)[1]
    gamma1[k] <- coef(fit)[2]
  }

  names(gamma0) <- names(gamma1) <- colnames(features)

  list(
    gamma0 = gamma0,
    gamma1 = gamma1
  )
}


## Reference profile representation

.make_reference_profile <- function(gamma0, gamma1) {
  .validate_named_numeric(gamma0, "gamma0")
  .validate_named_numeric(gamma1, "gamma1")

  ## require same name set (order can differ)
  .validate_same_names(gamma0, gamma1)

  ## align gamma1 to gamma0 order
  feat_names <- names(gamma0)
  gamma1 <- .reorder_named(gamma1, feat_names, name = "gamma1")

  ref <- list(gamma0 = gamma0, gamma1 = gamma1)
  class(ref) <- "reference_profile"
  ref
}

#' Create a reference age-feature profile
#'
#' Create a reference profile describing expected feature values as a linear
#' function of age: \eqn{E[X_k \mid age] = \gamma_{0k} + \gamma_{1k} \cdot age}.
#'
#' @param gamma0 Named numeric vector of intercepts per feature.
#' @param gamma1 Named numeric vector of slopes per feature.
#'
#' @return An object of class `reference_profile`.
#' @export
reference_profile <- function(gamma0, gamma1) {
  .make_reference_profile(gamma0 = gamma0, gamma1 = gamma1)
}

#' Fit per-feature linear models of feature values versus age to obtain a
#' `reference_profile`.
#'
#' @param features A matrix/data.frame of feature values (samples x features).
#' @param age Numeric vector of ages, length equal to nrow(features).
#'
#' @return An object of class `reference_profile`.
#'
#' @importFrom stats lm coef
#' @export
fit_reference_profile <- function(features, age) {
  trends <- .fit_reference_trends(features = features, age = age)
  .make_reference_profile(gamma0 = trends$gamma0, gamma1 = trends$gamma1)
}

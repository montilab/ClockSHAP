test_that("ClockSHAP contributions sum to deviation (core additivity)", {

  set.seed(1)

  # synthetic data
  features <- matrix(rnorm(20), nrow = 5, ncol = 4)
  colnames(features) <- paste0("F", 1:4)
  age <- c(40, 50, 60, 70, 80)

  # linear clock
  clock <- linear_clock(
    alpha = 10,
    beta  = setNames(runif(4), colnames(features)),
    mu    = setNames(rep(0, 4), colnames(features)),
    sigma = setNames(rep(1, 4), colnames(features))
  )

  # reference profile
  ref <- reference_profile(
    gamma0 = setNames(rep(0, 4), colnames(features)),
    gamma1 = setNames(rep(0.1, 4), colnames(features))
  )

  cs <- clockshap(features, age, clock, ref)

  expect_equal(
    rowSums(cs$phi),
    cs$deviation,
    tolerance = 1e-10
  )
})


test_that("predicted == expected + deviation", {

  set.seed(2)

  features <- matrix(rnorm(18), nrow = 6, ncol = 3)
  colnames(features) <- c("A", "B", "C")
  age <- seq(30, 55, length.out = 6)

  clock <- linear_clock(
    alpha = 3.25,
    beta  = c(A = 0.5, B = -0.25, C = 0.1),
    mu    = c(A = 0, B = 0, C = 0),
    sigma = c(A = 1, B = 1, C = 1)
  )

  ref <- reference_profile(
    gamma0 = c(A = 0.2, B = -0.1, C = 0.0),
    gamma1 = c(A = 0.01, B = 0.02, C = -0.005)
  )

  cs <- clockshap(features, age, clock, ref)

  expect_equal(cs$predicted, cs$expected + cs$deviation, tolerance = 1e-12)
})


test_that("ClockSHAP matches a manual computation with non-trivial mu/sigma", {

  set.seed(3)

  features <- matrix(rnorm(12), nrow = 4, ncol = 3)
  colnames(features) <- c("X", "Y", "Z")
  age <- c(40, 50, 60, 70)

  beta <- c(X = 0.9, Y = -0.4, Z = 0.2)
  mu <- c(X = 1.0, Y = -2.0, Z = 0.5)
  sigma <- c(X = 2.0, Y = 0.5, Z = 1.5)
  alpha <- -7

  clock <- linear_clock(alpha = alpha, beta = beta, mu = mu, sigma = sigma)

  gamma0 <- c(X = 0.5, Y = 1.5, Z = -1.0)
  gamma1 <- c(X = 0.01, Y = -0.02, Z = 0.005)
  ref <- reference_profile(gamma0 = gamma0, gamma1 = gamma1)

  cs <- clockshap(features, age, clock, ref)

  # Manual predicted
  z <- sweep(features, 2, mu, "-")
  z <- sweep(z, 2, sigma, "/")
  predicted_manual <- as.numeric(alpha + z %*% beta)

  # Manual expected
  expected_x <- sweep(matrix(age, nrow = length(age), ncol = length(gamma0)),
                      2, gamma1, "*")
  expected_x <- sweep(expected_x, 2, gamma0, "+")
  colnames(expected_x) <- names(gamma0)

  z_ref <- sweep(expected_x, 2, mu, "-")
  z_ref <- sweep(z_ref, 2, sigma, "/")
  expected_manual <- as.numeric(alpha + z_ref %*% beta)

  deviation_manual <- predicted_manual - expected_manual
  phi_manual <- sweep(z - z_ref, 2, beta, "*")

  expect_equal(cs$predicted, predicted_manual, tolerance = 1e-12)
  expect_equal(cs$expected, expected_manual, tolerance = 1e-12)
  expect_equal(cs$deviation, deviation_manual, tolerance = 1e-12)
  expect_equal(cs$phi, phi_manual, tolerance = 1e-12)
})


test_that("ClockSHAP is invariant to feature column order (name-based alignment)", {

  set.seed(4)

  features <- matrix(rnorm(20), nrow = 5, ncol = 4)
  colnames(features) <- c("F1", "F2", "F3", "F4")
  age <- c(40, 50, 60, 70, 80)

  clock <- linear_clock(
    alpha = 10,
    beta  = c(F1 = 0.2, F2 = -0.1, F3 = 0.05, F4 = 0.3),
    mu    = c(F1 = 1, F2 = 2, F3 = 3, F4 = 4),
    sigma = c(F1 = 2, F2 = 2, F3 = 2, F4 = 2)
  )

  ref <- reference_profile(
    gamma0 = c(F1 = 0, F2 = 0, F3 = 0, F4 = 0),
    gamma1 = c(F1 = 0.1, F2 = 0.2, F3 = 0.3, F4 = 0.4)
  )

  cs1 <- clockshap(features, age, clock, ref)

  perm <- c("F3", "F1", "F4", "F2")
  cs2 <- clockshap(features[, perm], age, clock, ref)

  # Compare by column names to avoid assumptions about internal ordering
  expect_equal(cs2$predicted, cs1$predicted, tolerance = 1e-12)
  expect_equal(cs2$expected,  cs1$expected,  tolerance = 1e-12)
  expect_equal(cs2$deviation, cs1$deviation, tolerance = 1e-12)
  expect_equal(cs2$phi[, colnames(cs1$phi)], cs1$phi, tolerance = 1e-12)
})


test_that("ClockSHAP gives identical results for matrix vs data.frame inputs", {

  set.seed(5)

  features <- matrix(rnorm(15), nrow = 5, ncol = 3)
  colnames(features) <- c("A", "B", "C")
  age <- c(25, 35, 45, 55, 65)

  clock <- linear_clock(
    alpha = 0,
    beta  = c(A = 1, B = 2, C = 3),
    mu    = c(A = 0, B = 0, C = 0),
    sigma = c(A = 1, B = 1, C = 1)
  )

  ref <- reference_profile(
    gamma0 = c(A = 0, B = 0, C = 0),
    gamma1 = c(A = 0, B = 0, C = 0)
  )

  cs_mat <- clockshap(features, age, clock, ref)
  cs_df  <- clockshap(as.data.frame(features), age, clock, ref)

  expect_equal(cs_df$predicted, cs_mat$predicted, tolerance = 1e-12)
  expect_equal(cs_df$expected,  cs_mat$expected,  tolerance = 1e-12)
  expect_equal(cs_df$deviation, cs_mat$deviation, tolerance = 1e-12)
  expect_equal(cs_df$phi,       cs_mat$phi,       tolerance = 1e-12)
})


test_that("ClockSHAP errors on mismatched feature names", {

  features <- matrix(rnorm(10), nrow = 5, ncol = 2)
  colnames(features) <- c("A", "B")
  age <- c(50, 55, 60, 65, 70)

  clock <- linear_clock(
    alpha = 0,
    beta  = c(A = 1, C = 2),
    mu    = c(A = 0, C = 0),
    sigma = c(A = 1, C = 1)
  )

  ref <- reference_profile(
    gamma0 = c(A = 0, C = 0),
    gamma1 = c(A = 0, C = 0)
  )

  expect_error(clockshap(features, age, clock, ref))
})


test_that("ClockSHAP errors when age length != nrow(features)", {

  features <- matrix(rnorm(12), nrow = 4, ncol = 3)
  colnames(features) <- c("A", "B", "C")

  clock <- linear_clock(
    alpha = 0,
    beta  = c(A = 1, B = 1, C = 1),
    mu    = c(A = 0, B = 0, C = 0),
    sigma = c(A = 1, B = 1, C = 1)
  )

  ref <- reference_profile(
    gamma0 = c(A = 0, B = 0, C = 0),
    gamma1 = c(A = 0, B = 0, C = 0)
  )

  expect_error(clockshap(features, age = c(10, 20, 30), clock, ref))
})


test_that("ClockSHAP errors on NA feature values (explicitly enforce complete cases)", {

  features <- matrix(rnorm(12), nrow = 4, ncol = 3)
  colnames(features) <- c("A", "B", "C")
  features[2, 1] <- NA

  age <- c(20, 30, 40, 50)

  clock <- linear_clock(
    alpha = 0,
    beta  = c(A = 1, B = 1, C = 1),
    mu    = c(A = 0, B = 0, C = 0),
    sigma = c(A = 1, B = 1, C = 1)
  )

  ref <- reference_profile(
    gamma0 = c(A = 0, B = 0, C = 0),
    gamma1 = c(A = 0, B = 0, C = 0)
  )

  expect_error(clockshap(features, age, clock, ref))
})


test_that("linear_clock errors when sigma has zeros", {

  expect_error(
    linear_clock(
      alpha = 0,
      beta  = c(A = 1),
      mu    = c(A = 0),
      sigma = c(A = 0)
    )
  )
})


test_that("ClockSHAP returns a well-formed object", {

  features <- matrix(0, nrow = 3, ncol = 2)
  colnames(features) <- c("X", "Y")
  age <- c(30, 40, 50)

  clock <- linear_clock(
    alpha = 0,
    beta  = c(X = 1, Y = 1),
    mu    = c(X = 0, Y = 0),
    sigma = c(X = 1, Y = 1)
  )

  ref <- reference_profile(
    gamma0 = c(X = 0, Y = 0),
    gamma1 = c(X = 0, Y = 0)
  )

  cs <- clockshap(features, age, clock, ref)

  expect_s3_class(cs, "clockshap")
  expect_true(all(c("phi", "deviation", "predicted", "expected", "age") %in% names(cs)))

  expect_true(is.matrix(cs$phi))
  expect_equal(nrow(cs$phi), nrow(features))
  expect_equal(ncol(cs$phi), ncol(features))
  expect_equal(colnames(cs$phi), colnames(features))
})


test_that("ClockSHAP returns zero deviation and zero contributions when features match the reference", {

  age <- c(40, 40)

  gamma0 <- c(X = 1,  Y = -2)
  gamma1 <- c(X = 0.1, Y = 0.2)

  features <- rbind(
    gamma0 + gamma1 * age[1],
    gamma0 + gamma1 * age[2]
  )
  colnames(features) <- names(gamma0)

  clock <- linear_clock(
    alpha = 0,
    beta  = c(X = 1, Y = 1),
    mu    = c(X = 0, Y = 0),
    sigma = c(X = 1, Y = 1)
  )

  ref <- reference_profile(gamma0 = gamma0, gamma1 = gamma1)

  cs <- clockshap(features, age, clock, ref)

  expect_equal(cs$deviation, c(0, 0), tolerance = 1e-12)
  expect_true(all(cs$phi == 0))
})


test_that("fit_reference_profile recovers coefficients for perfectly linear data", {

  age <- c(10, 20, 30, 40, 50)

  gamma0 <- c(A = 1.5, B = -3)
  gamma1 <- c(A = 0.2, B = 0.05)

  features <- rbind(
    gamma0 + gamma1 * age[1],
    gamma0 + gamma1 * age[2],
    gamma0 + gamma1 * age[3],
    gamma0 + gamma1 * age[4],
    gamma0 + gamma1 * age[5]
  )
  colnames(features) <- names(gamma0)

  ref_hat <- fit_reference_profile(as.data.frame(features), age)

  expect_s3_class(ref_hat, "reference_profile")
  expect_equal(ref_hat$gamma0, gamma0, tolerance = 1e-12)
  expect_equal(ref_hat$gamma1, gamma1, tolerance = 1e-12)
})


test_that("Accessors return the expected components", {

  features <- matrix(rnorm(12), nrow = 4, ncol = 3)
  colnames(features) <- c("A", "B", "C")
  age <- c(20, 30, 40, 50)

  clock <- linear_clock(
    alpha = 1,
    beta  = c(A = 0.1, B = 0.2, C = 0.3),
    mu    = c(A = 0, B = 0, C = 0),
    sigma = c(A = 1, B = 1, C = 1)
  )

  ref <- reference_profile(
    gamma0 = c(A = 0, B = 0, C = 0),
    gamma1 = c(A = 0, B = 0, C = 0)
  )

  cs <- clockshap(features, age, clock, ref)

  expect_equal(clockshap_phi(cs), cs$phi)
  expect_equal(clockshap_predicted(cs), cs$predicted)
  expect_equal(clockshap_expected(cs), cs$expected)
  expect_equal(clockshap_deviation(cs), cs$deviation)
})


test_that("plot_clockshap_waterfall returns a ggplot object", {

  skip_if_not_installed("ggplot2")

  set.seed(6)

  features <- matrix(rnorm(20), nrow = 5, ncol = 4)
  colnames(features) <- paste0("F", 1:4)
  rownames(features) <- paste0("S", 1:5)
  age <- c(40, 50, 60, 70, 80)

  clock <- linear_clock(
    alpha = 10,
    beta  = setNames(runif(4), colnames(features)),
    mu    = setNames(rep(0, 4), colnames(features)),
    sigma = setNames(rep(1, 4), colnames(features))
  )

  ref <- reference_profile(
    gamma0 = setNames(rep(0, 4), colnames(features)),
    gamma1 = setNames(rep(0.1, 4), colnames(features))
  )

  cs <- clockshap(features, age, clock, ref)

  p1 <- plot_clockshap_waterfall(cs, sample = 1, top_n = 3)
  p2 <- plot_clockshap_waterfall(cs, sample = "S1", top_n = 3)

  expect_s3_class(p1, "ggplot")
  expect_s3_class(p2, "ggplot")
})

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


test_that("ClockSHAP enforces strict alignment of names(age) and rownames(features)", {

  features <- matrix(rnorm(12), nrow = 4, ncol = 3)
  colnames(features) <- c("A", "B", "C")
  rownames(features) <- c("S1", "S2", "S3", "S4")

  age_ok <- c(S1 = 10, S2 = 20, S3 = 30, S4 = 40)
  age_bad_order <- c(S2 = 20, S1 = 10, S3 = 30, S4 = 40)
  age_bad_names <- c(S1 = 10, S2 = 20, S3 = 30, SX = 40)
  age_unnamed <- c(10, 20, 30, 40)

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

  expect_no_error(clockshap(features, age_ok, clock, ref))
  expect_error(clockshap(features, age_bad_order, clock, ref))
  expect_error(clockshap(features, age_bad_names, clock, ref))
  expect_error(clockshap(features, age_unnamed, clock, ref))
})


test_that("ClockSHAP errors if age is named but features rownames are missing", {

  features <- matrix(rnorm(12), nrow = 4, ncol = 3)
  colnames(features) <- c("A", "B", "C")
  rownames(features) <- NULL
  age <- c(S1 = 10, S2 = 20, S3 = 30, S4 = 40)

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


test_that("max_age marks out-of-range samples as NA in reference-based outputs", {

  features <- matrix(
    c(
      0.0,  1.0,
      0.5, -0.5,
      1.0,  0.0
    ),
    nrow = 3,
    byrow = TRUE
  )
  colnames(features) <- c("A", "B")
  age <- c(60, 75, 90)

  clock <- linear_clock(
    alpha = 5,
    beta  = c(A = 0.2, B = -0.1),
    mu    = c(A = 0, B = 0),
    sigma = c(A = 1, B = 1)
  )

  ref <- reference_profile(
    gamma0 = c(A = 0.0, B = 0.0),
    gamma1 = c(A = 0.01, B = 0.02)
  )

  cs_full <- clockshap(features, age, clock, ref)
  cs_cap  <- clockshap(features, age, clock, ref, max_age = 80)

  # Predicted is independent of reference age handling.
  expect_equal(cs_cap$predicted, cs_full$predicted, tolerance = 1e-12)

  # In-range rows stay finite and keep exact additivity.
  expect_false(anyNA(cs_cap$expected[1:2]))
  expect_false(anyNA(cs_cap$deviation[1:2]))
  expect_false(anyNA(cs_cap$phi[1:2, ]))
  expect_equal(rowSums(cs_cap$phi[1:2, , drop = FALSE]),
               cs_cap$deviation[1:2], tolerance = 1e-12)

  # Out-of-range rows are NA for reference-based outputs.
  expect_true(is.na(cs_cap$expected[3]))
  expect_true(is.na(cs_cap$deviation[3]))
  expect_true(all(is.na(cs_cap$phi[3, ])))
})


test_that("clockshap validates max_age", {

  features <- matrix(rnorm(8), nrow = 4, ncol = 2)
  colnames(features) <- c("A", "B")
  age <- c(30, 40, 50, 60)

  clock <- linear_clock(
    alpha = 0,
    beta  = c(A = 1, B = 1),
    mu    = c(A = 0, B = 0),
    sigma = c(A = 1, B = 1)
  )

  ref <- reference_profile(
    gamma0 = c(A = 0, B = 0),
    gamma1 = c(A = 0, B = 0)
  )

  expect_error(clockshap(features, age, clock, ref, max_age = NA_real_))
  expect_error(clockshap(features, age, clock, ref, max_age = c(70, 80)))
  expect_error(clockshap(features, age, clock, ref, max_age = "80"))
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


test_that("clockshap_example runs end to end", {

  data(clockshap_example, package = "ClockSHAP")

  expect_named(clockshap_example, c("features", "age", "clock"))
  expect_equal(dim(clockshap_example$features), c(120, 8))
  expect_identical(names(clockshap_example$age), rownames(clockshap_example$features))
  expect_s3_class(clockshap_example$clock, "linear_clock")

  ref <- fit_reference_profile(
    features = as.data.frame(clockshap_example$features),
    age = clockshap_example$age
  )
  cs <- clockshap(
    features = clockshap_example$features,
    age = clockshap_example$age,
    clock = clockshap_example$clock,
    reference = ref
  )

  expect_gt(cor(cs$predicted, clockshap_example$age), 0.8)
  expect_equal(unname(rowSums(cs$phi)), unname(cs$deviation), tolerance = 1e-10)
  expect_gt(cs$deviation[1], 0)
  expect_lt(cs$deviation[2], 0)
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
  age <- c(S1 = 40, S2 = 50, S3 = 60, S4 = 70, S5 = 80)

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
  p3 <- plot_clockshap_waterfall(
    cs,
    sample = 1,
    top_n = 3,
    x_axis_title = "Clock Age",
    y_axis_title = "Features",
    deviation_label = "Dev = {value} years",
    expected_age_label = "Expected\\n{value} years",
    predicted_age_label = "Predicted\\n{value} years",
    effect_label = "{value}"
  )

  expect_s3_class(p1, "ggplot")
  expect_s3_class(p2, "ggplot")
  expect_s3_class(p3, "ggplot")
  expect_no_warning(print(p1))
  expect_no_warning(print(p2))
})


test_that("waterfall effect labels match per-bar effects", {

  skip_if_not_installed("ggplot2")

  features <- matrix(
    c(
      1, 2, 3, 4,
      0, 0, 0, 0,
      0, 0, 0, 0
    ),
    nrow = 3,
    byrow = TRUE
  )
  colnames(features) <- paste0("F", 1:4)
  age <- c(40, 50, 60)

  clock <- linear_clock(
    alpha = 0,
    beta  = c(F1 = 1, F2 = 2, F3 = 3, F4 = 4),
    mu    = c(F1 = 0, F2 = 0, F3 = 0, F4 = 0),
    sigma = c(F1 = 1, F2 = 1, F3 = 1, F4 = 1)
  )

  ref <- reference_profile(
    gamma0 = c(F1 = 0, F2 = 0, F3 = 0, F4 = 0),
    gamma1 = c(F1 = 0, F2 = 0, F3 = 0, F4 = 0)
  )

  cs <- clockshap(features, age, clock, ref)

  p <- plot_clockshap_waterfall(
    cs,
    sample = 1,
    top_n = 3,
    show_age_labels = FALSE,
    show_delta_label = FALSE,
    effect_label = "phi={value}"
  )

  phi_all <- clockshap_phi(cs)[1, ]
  phi_all <- phi_all[order(-abs(phi_all))]
  top_eff <- utils::head(phi_all, 3)
  effects <- c(top_eff, Other = sum(phi_all) - sum(top_eff))
  expected_labels <- paste0("phi=", sprintf("%+.1f", effects))

  text_layer_idx <- which(
    vapply(p$layers, function(layer) inherits(layer$geom, "GeomText"), logical(1))
  )
  expect_equal(length(text_layer_idx), 1)

  built <- ggplot2::ggplot_build(p)
  expect_equal(as.character(built$data[[text_layer_idx]]$label), unname(expected_labels))
})


test_that("waterfall labels use one reconciled precision", {

  skip_if_not_installed("ggplot2")

  phi <- matrix(c(3.439, 1.234), nrow = 1,
                dimnames = list("S1", c("F1", "F2")))
  cs <- structure(
    list(
      phi = phi,
      deviation = 4.673,
      predicted = 57.034,
      expected = 52.361,
      age = 46
    ),
    class = "clockshap"
  )

  p <- plot_clockshap_waterfall(
    cs,
    sample = "S1",
    top_n = 2,
    effect_label_mode = "inside",
    deviation_label = "D:{value}",
    expected_age_label = "X:{value}",
    predicted_age_label = "P:{value}",
    effect_label = "E:{value}",
    digits = 2
  )

  built <- ggplot2::ggplot_build(p)
  labels <- unlist(lapply(
    built$data,
    function(layer) {
      if ("label" %in% names(layer)) as.character(layer$label) else character()
    }
  ))

  expect_true(all(c(
    "D:+4.67", "X:52.36", "P:57.03",
    "E:+3.44", "E:+1.23", "E:0.00"
  ) %in% labels))
  expect_error(
    plot_clockshap_waterfall(cs, sample = "S1", digits = 1.5),
    "integer between 0 and 9",
    fixed = TRUE
  )
  expect_warning(
    plot_clockshap_waterfall(cs, sample = "S1", top_label_digits = 2),
    "deprecated",
    fixed = TRUE
  )
  expect_error(
    suppressWarnings(plot_clockshap_waterfall(
      cs,
      sample = "S1",
      top_label_digits = 2,
      effect_label_digits = 1
    )),
    "must have the same value",
    fixed = TRUE
  )
})


test_that("reconciled waterfall values preserve displayed identities", {

  set.seed(2026)

  for (iteration in seq_len(20)) {
    phi <- matrix(rnorm(24, sd = 8), nrow = 2)
    colnames(phi) <- paste0("F", seq_len(ncol(phi)))
    rownames(phi) <- c("S1", "S2")
    deviation <- rowSums(phi)
    expected <- runif(2, 30, 80)
    cs <- structure(
      list(
        phi = phi,
        deviation = deviation,
        predicted = expected + deviation,
        expected = expected,
        age = c(45, 65)
      ),
      class = "clockshap"
    )

    for (digits in 0:2) {
      scale <- 10^digits
      for (top_n in c(1, 3, 7, 12)) {
        displayed <- clockshap_waterfall_values(
          cs, sample = "S1", top_n = top_n, digits = digits
        )
        repeated <- clockshap_waterfall_values(
          cs, sample = "S1", top_n = top_n, digits = digits
        )

        phi_sorted <- phi["S1", ][order(-abs(phi["S1", ]))]
        top <- head(phi_sorted, top_n)
        exact_contributions <- c(
          top,
          Other = sum(phi_sorted) - sum(top)
        )

        expect_equal(
          sum(round(displayed$contributions * scale)),
          round(displayed$deviation * scale)
        )
        expect_equal(
          round((displayed$expected + displayed$deviation) * scale),
          round(displayed$predicted * scale)
        )
        expect_lte(
          abs(displayed$predicted - cs$predicted[1]),
          1 / scale + 1e-12
        )
        expect_true(all(
          abs(
            round(displayed$contributions * scale) -
              round(exact_contributions * scale)
          ) <= 1
        ))
        nonzero_display <- displayed$contributions != 0
        expect_true(all(
          sign(displayed$contributions[nonzero_display]) ==
            sign(exact_contributions[nonzero_display])
        ))
        expect_identical(displayed, repeated)
      }
    }
  }
})


test_that("reconciliation handles near-zero deviations with cancellation", {

  phi <- matrix(
    c(100.49, -100.49, 0.49, 0.49),
    nrow = 1,
    dimnames = list("S1", c("A", "B", "C", "D"))
  )
  deviation <- rowSums(phi)
  cs <- structure(
    list(
      phi = phi,
      deviation = deviation,
      predicted = 50 + deviation,
      expected = 50,
      age = 50
    ),
    class = "clockshap"
  )

  displayed <- clockshap_waterfall_values(cs, "S1", top_n = 4, digits = 0)

  expect_equal(sum(displayed$contributions), displayed$deviation)
  expect_equal(displayed$expected + displayed$deviation, displayed$predicted)
  expect_equal(displayed$deviation, 1)
  expect_true(all(abs(displayed$contributions - round(c(phi[1, ], Other = 0))) <= 1))
})


test_that("waterfall label templates may omit the value placeholder", {

  skip_if_not_installed("ggplot2")

  phi <- matrix(c(1.2, -0.2), nrow = 1,
                dimnames = list("S1", c("F1", "F2")))
  cs <- structure(
    list(
      phi = phi,
      deviation = 1,
      predicted = 51,
      expected = 50,
      age = 50
    ),
    class = "clockshap"
  )

  p <- plot_clockshap_waterfall(
    cs,
    "S1",
    top_n = 2,
    effect_label_mode = "inside",
    deviation_label = "Age-matched deviation",
    expected_age_label = "Expected",
    predicted_age_label = "Predicted",
    effect_label = "Contribution"
  )
  built <- ggplot2::ggplot_build(p)
  labels <- unlist(lapply(
    built$data,
    function(layer) {
      if ("label" %in% names(layer)) as.character(layer$label) else character()
    }
  ))

  expect_true(all(c(
    "Age-matched deviation", "Expected", "Predicted", "Contribution"
  ) %in% labels))
})


test_that("waterfall auto mode moves small labels outside with connectors", {

  skip_if_not_installed("ggplot2")

  features <- matrix(c(10, 0.10, 0.05, 0.02), nrow = 1)
  colnames(features) <- c("F1", "F2", "F3", "F4")
  age <- c(50)

  clock <- linear_clock(
    alpha = 0,
    beta  = c(F1 = 1, F2 = 1, F3 = 1, F4 = 1),
    mu    = c(F1 = 0, F2 = 0, F3 = 0, F4 = 0),
    sigma = c(F1 = 1, F2 = 1, F3 = 1, F4 = 1)
  )
  ref <- reference_profile(
    gamma0 = c(F1 = 0, F2 = 0, F3 = 0, F4 = 0),
    gamma1 = c(F1 = 0, F2 = 0, F3 = 0, F4 = 0)
  )

  cs <- clockshap(features, age, clock, ref)

  p <- plot_clockshap_waterfall(
    cs,
    sample = 1,
    top_n = 4,
    show_age_labels = FALSE,
    show_delta_label = FALSE,
    effect_label = "phi={value}",
    effect_label_mode = "auto",
    effect_label_min_frac = 0.10,
    effect_label_outside_nudge_frac = 0.03,
    show_effect_label_connectors = TRUE
  )

  built <- ggplot2::ggplot_build(p)

  text_layer_idx <- which(
    vapply(p$layers, function(layer) inherits(layer$geom, "GeomText"), logical(1))
  )
  hjust_values <- unlist(lapply(text_layer_idx, function(i) built$data[[i]]$hjust))
  expect_true(any(hjust_values %in% c(0, 1)))

  has_connector_layer <- any(vapply(
    built$data,
    function(d) {
      is.data.frame(d) &&
        "colour" %in% names(d) &&
        nrow(d) > 0 &&
        all(d$colour == "grey40")
    },
    logical(1)
  ))
  expect_true(has_connector_layer)
})


test_that("waterfall hide_small mode suppresses small-bar effect labels", {

  skip_if_not_installed("ggplot2")

  features <- matrix(c(10, 0.10, 0.05, 0.02), nrow = 1)
  colnames(features) <- c("F1", "F2", "F3", "F4")
  age <- c(50)

  clock <- linear_clock(
    alpha = 0,
    beta  = c(F1 = 1, F2 = 1, F3 = 1, F4 = 1),
    mu    = c(F1 = 0, F2 = 0, F3 = 0, F4 = 0),
    sigma = c(F1 = 1, F2 = 1, F3 = 1, F4 = 1)
  )
  ref <- reference_profile(
    gamma0 = c(F1 = 0, F2 = 0, F3 = 0, F4 = 0),
    gamma1 = c(F1 = 0, F2 = 0, F3 = 0, F4 = 0)
  )

  cs <- clockshap(features, age, clock, ref)

  p <- plot_clockshap_waterfall(
    cs,
    sample = 1,
    top_n = 4,
    show_age_labels = FALSE,
    show_delta_label = FALSE,
    effect_label = "phi={value}",
    effect_label_mode = "hide_small",
    effect_label_min_frac = 0.10
  )

  built <- ggplot2::ggplot_build(p)
  text_layer_idx <- which(
    vapply(p$layers, function(layer) inherits(layer$geom, "GeomText"), logical(1))
  )
  labels <- unlist(lapply(text_layer_idx, function(i) as.character(built$data[[i]]$label)))
  expect_equal(unname(labels), "phi=+10.0")
})

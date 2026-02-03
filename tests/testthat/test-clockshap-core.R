test_that("ClockSHAP contributions sum to deviation", {

  set.seed(1)

  ## synthetic data
  features <- matrix(rnorm(20), nrow = 5, ncol = 4)
  colnames(features) <- paste0("F", 1:4)
  age <- c(40, 50, 60, 70, 80)

  ## linear clock
  clock <- linear_clock(
    alpha = 10,
    beta  = setNames(runif(4), colnames(features)),
    mu    = setNames(rep(0, 4), colnames(features)),
    sigma = setNames(rep(1, 4), colnames(features))
  )

  ## reference profile
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

  expect_error(
    clockshap(features, age, clock, ref)
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
  expect_named(cs, c("phi", "deviation", "predicted", "expected", "age"))
})

# test_coverage_ESCR.R
library(testthat)

# ── Helper: generate a standard fit object ──────────────────────────────────
make_fit <- function(seed = 123, L = 1000, alpha = 0.05,
                     optimize_type = "P", verbose = FALSE) {
  set.seed(seed)
  n <- 50
  x <- seq(-5, 5, length.out = n)
  X <- cbind(1, x, x^2)
  theta_true <- c(-6, -3, 0.25)
  Y <- X %*% theta_true + rnorm(n, 0, 0.2)
  compute_bscb_conjugate(
    X = X, Y = Y, alpha = alpha, a = -5, b = 5, L = L,
    theta_true = theta_true, optimize_type = optimize_type,
    verbose = verbose
  )
}

# ── Test 1: Return value is 0 or 1 ──────────────────────────────────────────
test_that("coverage_ESCR returns 0 or 1", {
  fit <- make_fit()
  result <- coverage_ESCR(fit, optimize_type = "P", verbose = FALSE)
  expect_true(result %in% c(0, 1))
})

# ── Test 2: Return type is numeric scalar ────────────────────────────────────
test_that("coverage_ESCR returns a scalar", {
  fit <- make_fit()
  result <- coverage_ESCR(fit, optimize_type = "P", verbose = FALSE)
  expect_length(result, 1)
  expect_true(is.numeric(result))
})

# ── Test 3: All three optimize_type options run without error ────────────────
test_that("all optimize_type options work", {
  fit <- make_fit()
  skip_if_not_installed("DEoptim")
  expect_no_error(coverage_ESCR(fit, optimize_type = "P", verbose = FALSE))
  expect_no_error(coverage_ESCR(fit, optimize_type = "G", verbose = FALSE))
  expect_no_error(coverage_ESCR(fit, optimize_type = "D", verbose = FALSE))
})

# ── Test 4: P and G give the same coverage flag ──────────────────────────────
test_that("P and G methods agree on coverage flag", {
  fit <- make_fit(L = 2000)
  result_P <- coverage_ESCR(fit, optimize_type = "P", verbose = FALSE)
  result_G <- coverage_ESCR(fit, optimize_type = "G", verbose = FALSE)
  expect_equal(result_P, result_G)
})

# ── Test 5: Invalid optimize_type is rejected ────────────────────────────────
test_that("invalid optimize_type is rejected", {
  fit <- make_fit()
  expect_error(
    coverage_ESCR(fit, optimize_type = "X", verbose = FALSE)
  )
})

# ── Test 6: fit without theta_true causes error ──────────────────────────────
test_that("fit without theta_true causes error", {
  set.seed(123)
  n <- 50
  x <- seq(-5, 5, length.out = n)
  X <- cbind(1, x, x^2)
  Y <- X %*% c(-6, -3, 0.25) + rnorm(n, 0, 0.2)

  fit_no_truth <- compute_bscb_conjugate(
    X = X, Y = Y, alpha = 0.05, a = -5, b = 5, L = 500,
    theta_true = NULL, verbose = FALSE
  )

  expect_error(
    coverage_ESCR(fit_no_truth, optimize_type = "P", verbose = FALSE)
  )
})

# ── Test 7: verbose = FALSE produces no output ───────────────────────────────
test_that("verbose = FALSE produces no output", {
  fit <- make_fit()
  expect_silent(
    coverage_ESCR(fit, optimize_type = "P", verbose = FALSE)
  )
})



# ── Test 8: Reproducibility with same fit object ─────────────────────────────
test_that("results are reproducible given the same fit object", {
  fit <- make_fit(seed = 42, L = 1000)
  # coverage_ESCR is deterministic given a fixed fit object
  r1 <- coverage_ESCR(fit, optimize_type = "P", verbose = FALSE)
  r2 <- coverage_ESCR(fit, optimize_type = "P", verbose = FALSE)
  expect_equal(r1, r2)
})

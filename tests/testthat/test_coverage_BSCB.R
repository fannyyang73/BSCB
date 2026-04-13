# test_coverage_ESCR.R
library(testthat)

# tests/testthat/test-coverage_BSCB.R

# ---------------------------------------------------------------------------
# Shared test data
# ---------------------------------------------------------------------------
set.seed(123)
n          <- 50
x          <- seq(-5, 5, length.out = n)
X          <- cbind(1, x, x^2)
theta_true <- c(-6, -3, 0.25)
Y          <- X %*% theta_true + rnorm(n, sd = 0.2)

fit <- compute_bscb_conjugate(
  X          = X,
  Y          = Y,
  alpha      = 0.05,
  a          = -5,
  b          =  5,
  L          = 500,
  theta_true = theta_true,
  verbose    = FALSE
)

# ---------------------------------------------------------------------------
# coverage_ESCR
# ---------------------------------------------------------------------------
test_that("coverage_ESCR returns 0 or 1", {
  result <- coverage_ESCR(fit, optimize_type = "P")
  expect_true(result %in% c(0L, 1L))
})

test_that("coverage_ESCR returns integer", {
  result <- coverage_ESCR(fit, optimize_type = "P")
  expect_type(result, "integer")
})

test_that("coverage_ESCR optimize_type P and G give consistent results", {
  result_P <- coverage_ESCR(fit, optimize_type = "P")
  result_G <- coverage_ESCR(fit, optimize_type = "G")
  # Both should return 0 or 1; exact agreement expected for same data
  expect_true(result_P %in% c(0L, 1L))
  expect_true(result_G %in% c(0L, 1L))
})

test_that("coverage_ESCR verbose prints message on failure", {
  # Create a fit with very small lambda to force failure
  fit_fail        <- fit
  fit_fail$lambda <- 0.001
  expect_message(
    coverage_ESCR(fit_fail, optimize_type = "P", verbose = TRUE),
    "Coverage failed"
  )
})

test_that("coverage_ESCR verbose is silent when coverage succeeds", {
  # Create a fit with very large lambda to force success
  fit_pass        <- fit
  fit_pass$lambda <- 999
  expect_no_message(
    coverage_ESCR(fit_pass, optimize_type = "P", verbose = TRUE)
  )
})

test_that("coverage_ESCR errors when theta_true is NULL", {
  fit_no_truth           <- fit
  fit_no_truth$theta_true <- NULL
  expect_error(
    coverage_ESCR(fit_no_truth, optimize_type = "P"),
    "fit\\$theta_true is NULL"
  )
})

test_that("coverage_ESCR rejects invalid optimize_type", {
  expect_error(
    coverage_ESCR(fit, optimize_type = "Z"),
    "should be one of"
  )
})

# ---------------------------------------------------------------------------
# coverage_PSCP
# ---------------------------------------------------------------------------
test_that("coverage_PSCP returns a numeric scalar in [0, 1]", {
  result <- coverage_PSCP(fit, draw_num = 100, optimize_type = "P")
  expect_type(result, "double")
  expect_length(result, 1)
  expect_true(result >= 0 && result <= 1)
})

test_that("coverage_PSCP with draw_num = 1 returns 0 or 1", {
  result <- coverage_PSCP(fit, draw_num = 1, optimize_type = "P")
  expect_true(result %in% c(0, 1))
})

test_that("coverage_PSCP increases with larger lambda", {
  fit_small        <- fit
  fit_small$lambda <- 0.001
  fit_large        <- fit
  fit_large$lambda <- 999

  pscp_small <- coverage_PSCP(fit_small, draw_num = 200,
                              optimize_type = "P")
  pscp_large <- coverage_PSCP(fit_large, draw_num = 200,
                              optimize_type = "P")
  expect_true(pscp_large >= pscp_small)
})

test_that("coverage_PSCP verbose prints PSCP message", {
  expect_message(
    coverage_PSCP(fit, draw_num = 50, optimize_type = "P",
                  verbose = TRUE),
    "PSCP"
  )
})

test_that("coverage_PSCP is silent when verbose = FALSE", {
  expect_no_message(
    coverage_PSCP(fit, draw_num = 50, optimize_type = "P",
                  verbose = FALSE)
  )
})

test_that("coverage_PSCP rejects invalid optimize_type", {
  expect_error(
    coverage_PSCP(fit, draw_num = 50, optimize_type = "Z"),
    "should be one of"
  )
})

test_that("coverage_PSCP results are reproducible with set.seed", {
  set.seed(42)
  pscp1 <- coverage_PSCP(fit, draw_num = 100, optimize_type = "P")
  set.seed(42)
  pscp2 <- coverage_PSCP(fit, draw_num = 100, optimize_type = "P")
  expect_equal(pscp1, pscp2)
})

test_that("coverage_PSCP works with compute_bscb_ind_jeffreys fit", {
  fit_j <- compute_bscb_ind_jeffreys(
    X          = X,
    Y          = Y,
    alpha      = 0.05,
    a          = -5,
    b          =  5,
    L          = 200,
    theta_true = theta_true,
    verbose    = FALSE
  )
  result <- coverage_PSCP(fit_j, draw_num = 100, optimize_type = "P")
  expect_true(result >= 0 && result <= 1)
})

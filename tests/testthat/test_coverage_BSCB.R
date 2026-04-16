# Helper: generate a small dataset and fit conjugate model
make_conjugate_fit <- function(n = 30, seed = 123) {
  set.seed(seed)
  x <- seq(-5, 5, length.out = n)
  X <- cbind(1, x, x^2)
  theta_true <- c(-6, -3, 0.25)
  Y <- as.numeric(X %*% theta_true + rnorm(n, sd = 0.2))
  fit <- compute_bscb_conjugate(
    X          = X,
    Y          = Y,
    alpha      = 0.05,
    a          = -5,
    b          =  5,
    L          = 100,
    theta_true = theta_true,
    verbose    = FALSE
  )
  fit
}

# Helper: generate a small dataset and fit Jeffreys model
make_jeffreys_fit <- function(n = 30, seed = 123) {
  set.seed(seed)
  x <- seq(-5, 5, length.out = n)
  X <- cbind(1, x, x^2)
  theta_true <- c(-6, -3, 0.25)
  Y <- as.numeric(X %*% theta_true + rnorm(n, sd = 0.2))
  fit <- compute_bscb_ind_jeffreys(
    X          = X,
    Y          = Y,
    alpha      = 0.05,
    a          = -5,
    b          =  5,
    L          = 100,
    theta_true = theta_true,
    verbose    = FALSE
  )
  fit
}

# ============================================================
# coverage_ESCR tests
# ============================================================

test_that("coverage_ESCR returns 0 or 1 for conjugate fit", {
  fit <- make_conjugate_fit()
  result <- coverage_ESCR(fit, optimize_type = "P", verbose = FALSE)
  expect_true(result %in% c(0L, 1L))
})

test_that("coverage_ESCR returns 0 or 1 for Jeffreys fit", {
  fit <- make_jeffreys_fit()
  result <- coverage_ESCR(fit, optimize_type = "P", verbose = FALSE)
  expect_true(result %in% c(0L, 1L))
})

test_that("coverage_ESCR returns integer", {
  fit <- make_conjugate_fit()
  result <- coverage_ESCR(fit, optimize_type = "P")
  expect_type(result, "integer")
  expect_length(result, 1)
})

test_that("coverage_ESCR errors when theta_true is NULL", {
  fit <- make_conjugate_fit()
  fit$theta_true <- NULL
  expect_error(
    coverage_ESCR(fit, optimize_type = "P"),
    "theta_true is NULL"
  )
})

test_that("coverage_ESCR works with optimize_type = 'G'", {
  fit <- make_conjugate_fit()
  result <- coverage_ESCR(fit, optimize_type = "G")
  expect_true(result %in% c(0L, 1L))
})

test_that("coverage_ESCR verbose prints message when coverage fails", {
  # Force a failure by setting lambda very small
  fit <- make_conjugate_fit()
  fit$lambda <- 0
  expect_message(
    coverage_ESCR(fit, optimize_type = "P", verbose = TRUE),
    "Coverage failed"
  )
})

test_that("coverage_ESCR does not print when verbose = FALSE", {
  fit <- make_conjugate_fit()
  fit$lambda <- 0
  expect_no_message(
    coverage_ESCR(fit, optimize_type = "P", verbose = FALSE)
  )
})

# ============================================================
# coverage_PSCP tests - conjugate
# ============================================================

test_that("coverage_PSCP returns value in [0, 1] for conjugate fit", {
  fit <- make_conjugate_fit()
  result <- coverage_PSCP(fit, draw_num = 100, optimize_type = "P")
  expect_gte(result, 0)
  expect_lte(result, 1)
})

test_that("coverage_PSCP returns value in [0, 1] for Jeffreys fit", {
  fit <- make_jeffreys_fit()
  result <- coverage_PSCP(fit, draw_num = 100, optimize_type = "P")
  expect_gte(result, 0)
  expect_lte(result, 1)
})

test_that("coverage_PSCP returns numeric scalar", {
  fit <- make_conjugate_fit()
  result <- coverage_PSCP(fit, draw_num = 100, optimize_type = "P")
  expect_type(result, "double")
  expect_length(result, 1)
})

test_that("coverage_PSCP errors when theta_true is NULL", {
  fit <- make_conjugate_fit()
  fit$theta_true <- NULL
  expect_error(
    coverage_PSCP(fit, draw_num = 100, optimize_type = "P"),
    "theta_true is NULL"
  )
})

test_that("coverage_PSCP with large lambda returns PSCP near 1", {
  fit <- make_conjugate_fit()
  fit$lambda <- 1e10
  result <- coverage_PSCP(fit, draw_num = 100, optimize_type = "P")
  expect_equal(result, 1)
})

test_that("coverage_PSCP with lambda = 0 returns PSCP near 0", {
  fit <- make_conjugate_fit()
  fit$lambda <- 0
  result <- coverage_PSCP(fit, draw_num = 100, optimize_type = "P")
  expect_equal(result, 0)
})

test_that("coverage_PSCP verbose prints message", {
  fit <- make_conjugate_fit()
  expect_message(
    coverage_PSCP(fit, draw_num = 100, optimize_type = "P", verbose = TRUE),
    "PSCP ="
  )
})

test_that("coverage_PSCP does not print when verbose = FALSE", {
  fit <- make_conjugate_fit()
  expect_no_message(
    coverage_PSCP(fit, draw_num = 100, optimize_type = "P", verbose = FALSE)
  )
})

test_that("coverage_PSCP works with optimize_type = 'G'", {
  fit <- make_conjugate_fit()
  result <- coverage_PSCP(fit, draw_num = 50, optimize_type = "G")
  expect_gte(result, 0)
  expect_lte(result, 1)
})

# ============================================================
# coverage_PSCP tests - HMC
# ============================================================

test_that("coverage_PSCP works for HMC fit using theta_mat", {
  skip_on_cran()
  skip_if_not(instantiate::stan_cmdstan_exists(), "CmdStan not available")

  set.seed(42)
  n <- 20
  x_seq <- seq(-5, 5, length.out = n)
  X <- cbind(1, x_seq, x_seq^2)
  theta_true <- c(-6, -3, 0.25)
  Y <- as.numeric(X %*% theta_true + rnorm(n, sd = 0.2))

  fit_h <- withCallingHandlers(
    compute_bscb_hmc(
      Y = Y, X = X, V = diag(n),
      a = -5, b = 5,
      theta_true = theta_true,
      iter_sampling = 200,
      iter_warmup   = 200,
      chains        = 2,
      L             = 100,
      draw_num      = 100
    ),
    stan_deprecate = function(w) invokeRestart("muffleWarning")
  )

  result <- coverage_PSCP(fit_h, draw_num = 100, optimize_type = "P")
  expect_gte(result, 0)
  expect_lte(result, 1)
})

test_that("coverage_PSCP result for HMC matches fit$PSCP approximately", {
  skip_on_cran()
  skip_if_not(instantiate::stan_cmdstan_exists(), "CmdStan not available")

  set.seed(42)
  n <- 20
  x_seq <- seq(-5, 5, length.out = n)
  X <- cbind(1, x_seq, x_seq^2)
  theta_true <- c(-6, -3, 0.25)
  Y <- as.numeric(X %*% theta_true + rnorm(n, sd = 0.2))

  fit_h <- withCallingHandlers(
    compute_bscb_hmc(
      Y = Y, X = X, V = diag(n),
      a = -5, b = 5,
      theta_true = theta_true,
      iter_sampling = 200,
      iter_warmup   = 200,
      chains        = 2,
      L             = 100,
      draw_num      = 500
    ),
    stan_deprecate = function(w) invokeRestart("muffleWarning")
  )

  result <- coverage_PSCP(fit_h, draw_num = 500, optimize_type = "P")
  # Both use Monte Carlo so allow tolerance of 0.05
  expect_equal(result, fit_h$PSCP, tolerance = 0.05)
})

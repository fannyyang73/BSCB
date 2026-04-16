
# ---------------------------------------------------------------------------
# Shared test data
# ---------------------------------------------------------------------------
set.seed(123)
n          <- 50
x          <- seq(-5, 5, length.out = n)
X          <- cbind(1, x, x^2)
theta_true <- c(-6, -3, 0.25)
Y          <- X %*% theta_true + rnorm(n, sd = 0.2)

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

# ---------------------------------------------------------------------------
# Output structure
# ---------------------------------------------------------------------------
test_that("returns a bscb_fit object with correct structure", {
  expect_s3_class(fit, "bscb_fit")
  expect_named(fit, c("lambda", "lower_bound", "upper_bound",
                      "theta_true", "order_form",
                      "mu_star",  "dof",
                      "scale_mat", "cov_theta",
                      "x_range", "call", "method",
                      "n", "p", "alpha",
                      "data", "lambda_samples", "params"))
})

test_that("method field is correct", {
  expect_equal(fit$method, "independent_jeffreys")
})

# ---------------------------------------------------------------------------
# Posterior parameters
# ---------------------------------------------------------------------------
test_that("mu_star has correct length", {
  expect_length(fit$mu_star, ncol(X))  # p + 1 = 3
})

test_that("cov_theta is a symmetric positive definite matrix", {
  expect_equal(nrow(fit$cov_theta), ncol(X))
  expect_equal(ncol(fit$cov_theta), ncol(X))
  expect_true(isSymmetric(fit$cov_theta))
  expect_true(all(eigen(fit$cov_theta)$values > 0))
})

test_that("dof equals n - p - 1", {
  expect_equal(fit$dof, n - 2 - 1)  # p = 2
})

# ---------------------------------------------------------------------------
# Critical constant lambda
# ---------------------------------------------------------------------------
test_that("lambda is a positive scalar", {
  expect_length(fit$lambda, 1)
  expect_true(fit$lambda > 0)
})

test_that("lambda_samples has length L", {
  expect_length(fit$lambda_samples, 100)
})

test_that("lambda equals the (1 - alpha) quantile of lambda_samples", {
  expect_equal(fit$lambda,
               as.numeric(quantile(fit$lambda_samples, probs = 0.95)),
               tolerance = 1e-10)
})

# ---------------------------------------------------------------------------
# Band functions
# ---------------------------------------------------------------------------
test_that("lower_bound and upper_bound are functions", {
  expect_true(is.function(fit$lower_bound))
  expect_true(is.function(fit$upper_bound))
})

test_that("lower_bound < upper_bound at all x values", {
  x_seq <- seq(-5, 5, length.out = 100)
  expect_true(all(fit$lower_bound(x_seq) < fit$upper_bound(x_seq)))
})

test_that("band functions return correct length for vector input", {
  x_seq <- seq(-5, 5, length.out = 100)
  expect_length(fit$lower_bound(x_seq), 100)
  expect_length(fit$upper_bound(x_seq), 100)
})

test_that("band functions return scalar for single x input", {
  expect_length(fit$lower_bound(0), 1)
  expect_length(fit$upper_bound(0), 1)
})

test_that("scalar and vectorized evaluation are consistent", {
  x_test <- c(-3, 0, 3)
  lower_vec    <- fit$lower_bound(x_test)
  lower_scalar <- sapply(x_test, fit$lower_bound)
  expect_equal(lower_vec, lower_scalar, tolerance = 1e-10)
})

# ---------------------------------------------------------------------------
# Input validation
# ---------------------------------------------------------------------------
test_that("error when Y length does not match nrow(X)", {
  expect_error(
    compute_bscb_ind_jeffreys(X = X, Y = Y[-1], L = 100, verbose = FALSE),
    "Length of Y must equal number of rows in X"
  )
})

test_that("error when AR_setting = 1 but rho is NULL", {
  expect_error(
    compute_bscb_ind_jeffreys(X = X, Y = Y, AR_setting = 1,
                              rho = NULL, L = 100, verbose = FALSE),
    "rho must be provided when AR_setting = 1"
  )
})

test_that("error for invalid AR_setting", {
  expect_error(
    compute_bscb_ind_jeffreys(X = X, Y = Y, AR_setting = 2,
                              L = 100, verbose = FALSE),
    "AR_setting must be 0 or 1"
  )
})

test_that("non-matrix X is coerced to matrix", {
  X_df <- as.data.frame(X)
  expect_no_error(
    compute_bscb_ind_jeffreys(X = X_df, Y = Y, L = 100, verbose = FALSE)
  )
})

# ---------------------------------------------------------------------------
# AR(1) errors
# ---------------------------------------------------------------------------
test_that("AR(1) setting runs without error", {
  expect_no_error(
    compute_bscb_ind_jeffreys(
      X          = X,
      Y          = Y,
      AR_setting = 1,
      rho        = 0.3,
      L          = 100,
      verbose    = FALSE
    )
  )
})

# ---------------------------------------------------------------------------
# x_range inference
# ---------------------------------------------------------------------------
test_that("x_range is inferred from X when a and b are NULL", {
  fit_inferred <- compute_bscb_ind_jeffreys(
    X = X, Y = Y, a = NULL, b = NULL, L = 100, verbose = FALSE
  )
  expect_equal(fit_inferred$x_range, c(min(X[, 2]), max(X[, 2])))
})

# ---------------------------------------------------------------------------
# optimize_type options
# ---------------------------------------------------------------------------
test_that("optimize_type = 'G' runs without error", {
  expect_no_error(
    compute_bscb_ind_jeffreys(
      X = X, Y = Y, L = 50,
      optimize_type = "G", verbose = FALSE
    )
  )
})

# ---------------------------------------------------------------------------
# scale_mat
# ---------------------------------------------------------------------------
test_that("scale_mat has correct dimension and relates to cov_theta", {
  expect_equal(dim(fit$scale_mat), c(ncol(X), ncol(X)))
  expect_equal(fit$cov_theta,
               (fit$dof / (fit$dof - 2)) * fit$scale_mat,
               tolerance = 1e-10)
})

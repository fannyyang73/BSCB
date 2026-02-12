# test_compute_bscb_conjugate.R
library(testthat)

# Test 1: Basic functionality with quadratic model
test_that("compute_bscb_conjugate works with quadratic model", {
  set.seed(123)
  n <- 50
  p <- 2
  x <- seq(-5, 5, length.out = n)
  X <- cbind(1, x, x^2)
  theta_true <- c(-6, -3, 0.25)
  Y <- X %*% theta_true + rnorm(n, 0, 0.2)

  fit <- compute_bscb_conjugate(
    X = X, Y = Y, alpha = 0.05, a = -5, b = 5, L = 1000,
    AR_setting = 0, verbose = FALSE
  )

  # Check structure
  expect_s3_class(fit, "bscb_fit")
  expect_true(is.numeric(fit$lambda))
  expect_true(is.function(fit$lower_bound))
  expect_true(is.function(fit$upper_bound))
  expect_equal(length(fit$mu_star), p + 1)
  expect_equal(dim(fit$cov_theta), c(p + 1, p + 1))

  # Check lambda is positive
  expect_true(fit$lambda > 0)

  # Check bounds work at a point
  x_test <- 0
  lower <- fit$lower_bound(x_test)
  upper <- fit$upper_bound(x_test)
  expect_true(is.numeric(lower))
  expect_true(is.numeric(upper))
  expect_true(upper > lower)
})

# Test 2: Vectorized bound functions
test_that("bound functions are vectorized", {
  set.seed(123)
  n <- 50
  x <- seq(-5, 5, length.out = n)
  X <- cbind(1, x, x^2)
  theta_true <- c(-6, -3, 0.25)
  Y <- X %*% theta_true + rnorm(n, 0, 0.2)

  fit <- compute_bscb_conjugate(
    X = X, Y = Y, alpha = 0.05, a = -5, b = 5, L = 1000,
    verbose = FALSE
  )

  # Test vectorization
  x_vec <- seq(-5, 5, length.out = 10)
  lower_vec <- fit$lower_bound(x_vec)
  upper_vec <- fit$upper_bound(x_vec)

  expect_equal(length(lower_vec), length(x_vec))
  expect_equal(length(upper_vec), length(x_vec))
  expect_true(all(upper_vec > lower_vec))
})

# Test 3: Different alpha levels
test_that("different alpha levels produce different lambdas", {
  set.seed(123)
  n <- 50
  x <- seq(-5, 5, length.out = n)
  X <- cbind(1, x, x^2)
  theta_true <- c(-6, -3, 0.25)
  Y <- X %*% theta_true + rnorm(n, 0, 0.2)

  fit_90 <- compute_bscb_conjugate(
    X = X, Y = Y, alpha = 0.10, a = -5, b = 5, L = 1000,
    verbose = FALSE
  )

  fit_95 <- compute_bscb_conjugate(
    X = X, Y = Y, alpha = 0.05, a = -5, b = 5, L = 1000,
    verbose = FALSE
  )

  fit_99 <- compute_bscb_conjugate(
    X = X, Y = Y, alpha = 0.01, a = -5, b = 5, L = 1000,
    verbose = FALSE
  )

  # Higher confidence should have larger lambda
  expect_true(fit_90$lambda < fit_95$lambda)
  expect_true(fit_95$lambda < fit_99$lambda)
})

# Test 4: AR errors
test_that("AR error structure works", {
  set.seed(456)
  n <- 50
  x <- seq(-5, 5, length.out = n)
  X <- cbind(1, x, x^2)
  theta_true <- c(-6, -3, 0.25)
  rho <- 0.5

  # Generate AR errors
  epsilon_ar <- numeric(n)
  epsilon_ar[1] <- rnorm(1, 0, 0.2)
  for (i in 2:n) {
    epsilon_ar[i] <- rho * epsilon_ar[i-1] + rnorm(1, 0, 0.2 * sqrt(1 - rho^2))
  }
  Y_ar <- X %*% theta_true + epsilon_ar

  fit_ar <- compute_bscb_conjugate(
    X = X, Y = Y_ar, alpha = 0.05, a = -5, b = 5, L = 1000,
    AR_setting = 1, rho = 0.5, verbose = FALSE
  )

  expect_s3_class(fit_ar, "bscb_fit")
  expect_equal(fit_ar$params$AR_setting, 1)
  expect_equal(fit_ar$params$rho, 0.5)
})

# Test 5: AR setting requires rho
test_that("AR setting without rho throws error", {
  set.seed(123)
  n <- 50
  x <- seq(-5, 5, length.out = n)
  X <- cbind(1, x, x^2)
  Y <- rnorm(n)

  expect_error(
    compute_bscb_conjugate(
      X = X, Y = Y, alpha = 0.05, a = -5, b = 5, L = 1000,
      AR_setting = 1, rho = NULL, verbose = FALSE
    ),
    "rho must be provided when AR_setting = 1"
  )
})

# Test 6: Input validation
test_that("input validation works", {
  set.seed(123)
  n <- 50
  x <- seq(-5, 5, length.out = n)
  X <- cbind(1, x, x^2)
  Y <- rnorm(n)

  # Y length mismatch
  expect_error(
    compute_bscb_conjugate(
      X = X, Y = rnorm(n + 5), alpha = 0.05, a = -5, b = 5, L = 1000,
      verbose = FALSE
    ),
    "Length of Y must equal number of rows in X"
  )

  # Invalid AR_setting
  expect_error(
    compute_bscb_conjugate(
      X = X, Y = Y, alpha = 0.05, a = -5, b = 5, L = 1000,
      AR_setting = 2, verbose = FALSE
    ),
    "AR_setting must be 0 or 1"
  )
})

# Test 7: Automatic range detection
test_that("automatic x range detection works", {
  set.seed(123)
  n <- 50
  x <- seq(-3, 7, length.out = n)
  X <- cbind(1, x, x^2)
  Y <- rnorm(n)

  fit <- compute_bscb_conjugate(
    X = X, Y = Y, alpha = 0.05, L = 1000,
    verbose = FALSE
  )

  expect_equal(fit$x_range[1], min(x))
  expect_equal(fit$x_range[2], max(x))
})

# Test 8: Different polynomial degrees
test_that("works with different polynomial degrees", {
  set.seed(123)
  n <- 100

  # Linear (p=1)
  x <- seq(-5, 5, length.out = n)
  X_linear <- cbind(1, x)
  theta_linear <- c(2, 0.5)
  Y_linear <- X_linear %*% theta_linear + rnorm(n, 0, 0.2)

  fit_linear <- compute_bscb_conjugate(
    X = X_linear, Y = Y_linear, alpha = 0.05, a = -5, b = 5, L = 1000,
    verbose = FALSE
  )

  expect_equal(fit_linear$p, 1)
  expect_equal(length(fit_linear$mu_star), 2)

  # Cubic (p=3)
  X_cubic <- cbind(1, x, x^2, x^3)
  theta_cubic <- c(1, 0, -1, 0.2)
  Y_cubic <- X_cubic %*% theta_cubic + rnorm(n, 0, 0.5)

  fit_cubic <- compute_bscb_conjugate(
    X = X_cubic, Y = Y_cubic, alpha = 0.05, a = -5, b = 5, L = 1000,
    verbose = FALSE
  )

  expect_equal(fit_cubic$p, 3)
  expect_equal(length(fit_cubic$mu_star), 4)
})

# Test 9: Different optimization methods
test_that("both optimization methods work", {
  set.seed(123)
  n <- 50
  x <- seq(-5, 5, length.out = n)
  X <- cbind(1, x, x^2)
  theta_true <- c(-6, -3, 0.25)
  Y <- X %*% theta_true + rnorm(n, 0, 0.2)

  # Skip DEoptim test if package not available
  skip_if_not_installed("DEoptim")

  fit_D <- compute_bscb_conjugate(
    X = X, Y = Y, alpha = 0.05, a = -5, b = 5, L = 500,
    optimize_type = "D", verbose = FALSE
  )

  fit_G <- compute_bscb_conjugate(
    X = X, Y = Y, alpha = 0.05, a = -5, b = 5, L = 500,
    optimize_type = "G", verbose = FALSE
  )

  expect_s3_class(fit_D, "bscb_fit")
  expect_s3_class(fit_G, "bscb_fit")

  # Both should produce reasonable lambdas
  expect_true(fit_D$lambda > 0)
  expect_true(fit_G$lambda > 0)
})

# Test 10: Return object structure
test_that("return object has correct structure", {
  set.seed(123)
  n <- 50
  x <- seq(-5, 5, length.out = n)
  X <- cbind(1, x, x^2)
  Y <- rnorm(n)

  fit <- compute_bscb_conjugate(
    X = X, Y = Y, alpha = 0.05, a = -5, b = 5, L = 1000,
    verbose = FALSE
  )

  # Check all required fields
  expect_true("lambda" %in% names(fit))
  expect_true("lower_bound" %in% names(fit))
  expect_true("upper_bound" %in% names(fit))
  expect_true("mu_star" %in% names(fit))
  expect_true("cov_theta" %in% names(fit))
  expect_true("dof" %in% names(fit))
  expect_true("x_range" %in% names(fit))
  expect_true("order_form" %in% names(fit))
  expect_true("theta_true" %in% names(fit))
  expect_true("method" %in% names(fit))
  expect_true("data" %in% names(fit))
  expect_true("params" %in% names(fit))
  expect_true("lambda_samples" %in% names(fit))

  # Check metadata
  expect_equal(fit$method, "conjugate")
  expect_equal(fit$n, n)
  expect_equal(fit$p, 2)
  expect_equal(fit$alpha, 0.05)
})

# Test 11: Lambda samples
test_that("lambda samples are stored correctly", {
  set.seed(123)
  n <- 50
  x <- seq(-5, 5, length.out = n)
  X <- cbind(1, x, x^2)
  Y <- rnorm(n)
  L <- 1000

  fit <- compute_bscb_conjugate(
    X = X, Y = Y, alpha = 0.05, a = -5, b = 5, L = L,
    verbose = FALSE
  )

  expect_equal(length(fit$lambda_samples), L)
  expect_true(all(fit$lambda_samples > 0))

  # Lambda should be the (1-alpha) quantile
  manual_lambda <- quantile(fit$lambda_samples, probs = 0.95)
  expect_equal(fit$lambda, as.numeric(manual_lambda))
})

# Test 12: Prior settings
test_that("different prior settings work", {
  set.seed(123)
  n <- 50
  x <- seq(-5, 5, length.out = n)
  X <- cbind(1, x, x^2)
  Y <- rnorm(n)

  # P_setting = 1
  fit_P1 <- compute_bscb_conjugate(
    X = X, Y = Y, alpha = 0.05, a = -5, b = 5, L = 1000,
    P_setting = 1, c_value = 0.001, verbose = FALSE
  )

  # P_setting = 2
  fit_P2 <- compute_bscb_conjugate(
    X = X, Y = Y, alpha = 0.05, a = -5, b = 5, L = 1000,
    P_setting = 2, verbose = FALSE
  )

  expect_equal(fit_P1$params$P_setting, 1)
  expect_equal(fit_P2$params$P_setting, 2)
  expect_s3_class(fit_P1, "bscb_fit")
  expect_s3_class(fit_P2, "bscb_fit")
})

# Test 13: Coverage simulation (integration test)
test_that("empirical coverage rate is reasonable", {
  skip_on_cran()  # This test is slow

  set.seed(123)
  n <- 50
  x <- seq(-5, 5, length.out = n)
  X <- cbind(1, x, x^2)
  theta_true <- c(-6, -3, 0.25)

  n_sim <- 50  # Small number for testing
  coverage_count <- 0

  for (i in 1:n_sim) {
    set.seed(1000 + i)
    Y <- X %*% theta_true + rnorm(n, 0, 0.2)

    fit <- compute_bscb_conjugate(
      X = X, Y = Y, alpha = 0.05, a = -5, b = 5, L = 5000,
      verbose = FALSE
    )

    # Check coverage at multiple points
    x_check <- seq(-5, 5, length.out = 20)
    X_check <- cbind(1, x_check, x_check^2)
    y_true <- X_check %*% theta_true
    lower <- fit$lower_bound(x_check)
    upper <- fit$upper_bound(x_check)

    covered <- all(y_true >= lower & y_true <= upper)
    coverage_count <- coverage_count + covered
  }

  coverage_rate <- coverage_count / n_sim

  # Should be approximately 0.95, allow wide tolerance for small n_sim
  expect_true(coverage_rate >= 0.80 && coverage_rate <= 1.0)
})

# Test 14: Consistency of results with same seed
test_that("results are reproducible with same seed", {
  n <- 50
  x <- seq(-5, 5, length.out = n)
  X <- cbind(1, x, x^2)
  Y <- rnorm(n)

  set.seed(123)
  fit1 <- compute_bscb_conjugate(
    X = X, Y = Y, alpha = 0.05, a = -5, b = 5, L = 1000,
    verbose = FALSE
  )

  set.seed(123)
  fit2 <- compute_bscb_conjugate(
    X = X, Y = Y, alpha = 0.05, a = -5, b = 5, L = 1000,
    verbose = FALSE
  )

  expect_equal(fit1$lambda, fit2$lambda)
  expect_equal(fit1$mu_star, fit2$mu_star)
  expect_equal(fit1$lambda_samples, fit2$lambda_samples)
})

# Test 15: Verbose output
test_that("verbose mode produces messages", {
  set.seed(123)
  n <- 50
  x <- seq(-5, 5, length.out = n)
  X <- cbind(1, x, x^2)
  Y <- rnorm(n)

  expect_message(
    compute_bscb_conjugate(
      X = X, Y = Y, alpha = 0.05, a = -5, b = 5, L = 100,
      verbose = TRUE
    ),
    "Computing lambda"
  )

  expect_message(
    compute_bscb_conjugate(
      X = X, Y = Y, alpha = 0.05, a = -5, b = 5, L = 100,
      verbose = TRUE
    ),
    "lambda ="
  )
})

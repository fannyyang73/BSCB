# Helper function to generate test data
generate_test_data <- function(n = 50, p = 2, seed = 123) {
  set.seed(seed)

  x <- seq(-5, 5, length.out = n)
  X <- cbind(1, x, x^2)
  theta_true <- c(-6, -3, 0.25)
  e_sd <- 0.2

  Y <- X %*% theta_true + rnorm(n, sd = e_sd)

  return(list(X = X, Y = Y, theta_true = theta_true, n = n, p = p))
}

test_that("compute_bscb_conjugate produces valid output", {
  # Generate test data
  test_data <- generate_test_data(n = 50, p = 2)

  fit <- compute_bscb_conjugate(
    X = test_data$X,
    Y = test_data$Y,
    alpha = 0.05,
    a = -5,
    b = 5,
    L = 100,  # 小数值用于快速测试
    AR_setting = 0,
    optimize_type = "G",
    verbose = FALSE
  )

  # Test structure
  expect_s3_class(fit, "bscb_fit")
  expect_true(!is.null(fit$lambda))
  expect_true(!is.null(fit$mu_star))
  expect_true(!is.null(fit$lower_bound))
  expect_true(!is.null(fit$upper_bound))

  # Test lambda is positive
  expect_true(fit$lambda > 0)

  # Test mu_star length
  expect_equal(length(fit$mu_star), test_data$p + 1)

  # Test bound functions
  x_new <- 0.5
  lower <- fit$lower_bound(x_new)
  upper <- fit$upper_bound(x_new)

  expect_true(is.numeric(lower))
  expect_true(is.numeric(upper))
  expect_true(lower < upper)

  # Test vectorization
  x_seq <- seq(-5, 5, length.out = 10)
  lower_vec <- fit$lower_bound(x_seq)
  upper_vec <- fit$upper_bound(x_seq)

  expect_equal(length(lower_vec), 10)
  expect_equal(length(upper_vec), 10)
  expect_true(all(lower_vec < upper_vec))
})

test_that("compute_bscb_conjugate handles AR errors", {
  test_data <- generate_test_data(n = 50, p = 2)

  fit_ar <- compute_bscb_conjugate(
    X = test_data$X,
    Y = test_data$Y,
    alpha = 0.05,
    a = -5,
    b = 5,
    L = 100,
    AR_setting = 1,
    rho = 0.5,
    optimize_type = "G",
    verbose = FALSE
  )

  expect_s3_class(fit_ar, "bscb_fit")
  expect_true(!is.null(fit_ar$lambda))
})

test_that("compute_bscb_conjugate errors without rho for AR", {
  test_data <- generate_test_data(n = 50, p = 2)

  expect_error(
    compute_bscb_conjugate(
      X = test_data$X,
      Y = test_data$Y,
      AR_setting = 1,
      rho = NULL
    ),
    "rho must be provided"
  )
})

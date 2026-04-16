# Helper: generate a small dataset for testing
make_test_data <- function(n = 20, p = 2, seed = 123) {
  set.seed(seed)
  x_seq <- seq(-5, 5, length.out = n)
  if (p == 2) {
    X          <- cbind(1, x_seq, x_seq^2)
    theta_true <- c(-6, -3, 0.25)
  } else {
    X          <- cbind(1, x_seq, x_seq^2, x_seq^3)
    theta_true <- c(1, 2, -1, 0.5)
  }
  Y <- as.numeric(X %*% theta_true + rnorm(n, sd = 0.2))
  list(X = X, Y = Y, theta_true = theta_true, n = n, p = p)
}

# Shared HMC settings: minimal iterations to keep tests fast
HMC_fast <- list(
  iter_sampling = 200,
  iter_warmup   = 200,
  chains        = 2,
  L             = 100,
  draw_num      = 100
)

# ============================================================
# 1. Return structure
# ============================================================
test_that("compute_bscb_hmc returns a bscb_fit object with correct fields (p=2)", {
  skip_on_cran()
  skip_if_not(instantiate::stan_cmdstan_exists(), "CmdStan not available")

  d   <- make_test_data(p = 2)
  fit <- do.call(compute_bscb_hmc, c(
    list(Y = d$Y, X = d$X, V = diag(d$n),
         a = -5, b = 5, theta_true = d$theta_true),
    HMC_fast
  ))

  expect_s3_class(fit, "bscb_fit")
  expect_named(fit, c(
    "lambda", "lower_bound", "upper_bound", "theta_true", "order_form",
    "PSCP", "mu_star", "cov_theta", "theta_mat", "x_range",
    "call", "method", "n", "p", "alpha", "data",
    "lambda_samples", "params"
  ), ignore.order = TRUE)
})

# ============================================================
# 2. Metadata correctness
# ============================================================
test_that("compute_bscb_hmc stores correct metadata", {
  skip_on_cran()
  skip_if_not(instantiate::stan_cmdstan_exists(), "CmdStan not available")

  d   <- make_test_data(p = 2)
  fit <- do.call(compute_bscb_hmc, c(
    list(Y = d$Y, X = d$X, V = diag(d$n),
         a = -5, b = 5, theta_true = d$theta_true, alpha = 0.05),
    HMC_fast
  ))

  expect_equal(fit$method, "HMC")
  expect_equal(fit$n, d$n)
  expect_equal(fit$p, d$p)
  expect_equal(fit$alpha, 0.05)
  expect_equal(fit$x_range, c(-5, 5))
  expect_equal(fit$theta_true, d$theta_true)
})

# ============================================================
# 3. lambda is a positive scalar
# ============================================================
test_that("lambda is a positive scalar", {
  skip_on_cran()
  skip_if_not(instantiate::stan_cmdstan_exists(), "CmdStan not available")

  d   <- make_test_data(p = 2)
  fit <- do.call(compute_bscb_hmc, c(
    list(Y = d$Y, X = d$X, V = diag(d$n), a = -5, b = 5),
    HMC_fast
  ))

  expect_true(is.numeric(fit$lambda))
  expect_length(fit$lambda, 1)
  expect_gt(fit$lambda, 0)
})

# ============================================================
# 4. lambda_samples has correct length
# ============================================================
test_that("lambda_samples has length L", {
  skip_on_cran()
  skip_if_not(instantiate::stan_cmdstan_exists(), "CmdStan not available")

  d   <- make_test_data(p = 2)
  fit <- do.call(compute_bscb_hmc, c(
    list(Y = d$Y, X = d$X, V = diag(d$n), a = -5, b = 5),
    HMC_fast
  ))

  expect_length(fit$lambda_samples, HMC_fast$L)
  expect_true(all(fit$lambda_samples > 0))
})

# ============================================================
# 5. Posterior dimensions
# ============================================================
test_that("mu_star and cov_theta have correct dimensions (p=2)", {
  skip_on_cran()
  skip_if_not(instantiate::stan_cmdstan_exists(), "CmdStan not available")

  d   <- make_test_data(p = 2)
  fit <- do.call(compute_bscb_hmc, c(
    list(Y = d$Y, X = d$X, V = diag(d$n), a = -5, b = 5),
    HMC_fast
  ))

  expect_length(fit$mu_star, d$p + 1)
  expect_equal(dim(fit$cov_theta), c(d$p + 1, d$p + 1))
  expect_true(isSymmetric(fit$cov_theta))
})

# ============================================================
# 6. theta_mat dimensions
# ============================================================
test_that("theta_mat has correct dimensions", {
  skip_on_cran()
  skip_if_not(instantiate::stan_cmdstan_exists(), "CmdStan not available")

  d   <- make_test_data(p = 2)
  fit <- do.call(compute_bscb_hmc, c(
    list(Y = d$Y, X = d$X, V = diag(d$n), a = -5, b = 5),
    HMC_fast
  ))

  expected_rows <- HMC_fast$iter_sampling * HMC_fast$chains
  expect_equal(ncol(fit$theta_mat), d$p + 1)
  expect_equal(nrow(fit$theta_mat), expected_rows)
})

# ============================================================
# 7. lower_bound and upper_bound are functions
# ============================================================
test_that("lower_bound and upper_bound are callable functions", {
  skip_on_cran()
  skip_if_not(instantiate::stan_cmdstan_exists(), "CmdStan not available")

  d   <- make_test_data(p = 2)
  fit <- do.call(compute_bscb_hmc, c(
    list(Y = d$Y, X = d$X, V = diag(d$n), a = -5, b = 5),
    HMC_fast
  ))

  expect_true(is.function(fit$lower_bound))
  expect_true(is.function(fit$upper_bound))
})

# ============================================================
# 8. Band ordering: lower < upper at every point
# ============================================================
test_that("lower_bound < upper_bound across [a, b]", {
  skip_on_cran()
  skip_if_not(instantiate::stan_cmdstan_exists(), "CmdStan not available")

  d   <- make_test_data(p = 2)
  fit <- do.call(compute_bscb_hmc, c(
    list(Y = d$Y, X = d$X, V = diag(d$n), a = -5, b = 5),
    HMC_fast
  ))

  x_grid <- seq(-5, 5, length.out = 100)
  expect_true(all(fit$lower_bound(x_grid) < fit$upper_bound(x_grid)))
})

# ============================================================
# 9. Vectorised band evaluation
# ============================================================
test_that("lower_bound and upper_bound return vectors of correct length", {
  skip_on_cran()
  skip_if_not(instantiate::stan_cmdstan_exists(), "CmdStan not available")

  d   <- make_test_data(p = 2)
  fit <- do.call(compute_bscb_hmc, c(
    list(Y = d$Y, X = d$X, V = diag(d$n), a = -5, b = 5),
    HMC_fast
  ))

  x_grid <- seq(-5, 5, length.out = 50)
  expect_length(fit$lower_bound(x_grid), 50)
  expect_length(fit$upper_bound(x_grid), 50)
})

# ============================================================
# 10. PSCP is in [0, 1]
# ============================================================
test_that("PSCP is a probability in [0, 1]", {
  skip_on_cran()
  skip_if_not(instantiate::stan_cmdstan_exists(), "CmdStan not available")

  d   <- make_test_data(p = 2)
  fit <- do.call(compute_bscb_hmc, c(
    list(Y = d$Y, X = d$X, V = diag(d$n), a = -5, b = 5),
    HMC_fast
  ))

  expect_gte(fit$PSCP, 0)
  expect_lte(fit$PSCP, 1)
})

# ============================================================
# 11. p = 3 (cubic) works correctly
# ============================================================
test_that("compute_bscb_hmc works for p=3 (cubic)", {
  skip_on_cran()
  skip_if_not(instantiate::stan_cmdstan_exists(), "CmdStan not available")

  d   <- make_test_data(p = 3)
  fit <- do.call(compute_bscb_hmc, c(
    list(Y = d$Y, X = d$X, V = diag(d$n), a = -5, b = 5),
    HMC_fast
  ))

  expect_s3_class(fit, "bscb_fit")
  expect_equal(fit$p, 3)
  expect_length(fit$mu_star, 4)
  expect_equal(dim(fit$cov_theta), c(4, 4))
})

# ============================================================
# 12. prior_type = "normal_normal" works
# ============================================================
test_that("compute_bscb_hmc works with normal_normal prior", {
  skip_on_cran()
  skip_if_not(instantiate::stan_cmdstan_exists(), "CmdStan not available")

  d   <- make_test_data(p = 2)
  fit <- do.call(compute_bscb_hmc, c(
    list(Y = d$Y, X = d$X, V = diag(d$n),
         a = -5, b = 5, prior_type = "normal_normal"),
    HMC_fast
  ))

  expect_s3_class(fit, "bscb_fit")
  expect_equal(fit$params$prior_type, "normal_normal")
  expect_gt(fit$lambda, 0)
})

# ============================================================
# 13. optimize_type = "G" works
# ============================================================
test_that("compute_bscb_hmc works with optimize_type = 'G'", {
  skip_on_cran()
  skip_if_not(instantiate::stan_cmdstan_exists(), "CmdStan not available")

  d   <- make_test_data(p = 2)
  fit <- do.call(compute_bscb_hmc, c(
    list(Y = d$Y, X = d$X, V = diag(d$n),
         a = -5, b = 5, optimize_type = "G"),
    HMC_fast
  ))

  expect_s3_class(fit, "bscb_fit")
  expect_gt(fit$lambda, 0)
})

# ============================================================
# 14. Input validation: mismatched Y and X
# ============================================================
test_that("compute_bscb_hmc errors on mismatched Y and X", {
  d <- make_test_data(p = 2)
  expect_error(
    compute_bscb_hmc(Y = d$Y[-1], X = d$X, V = diag(d$n), a = -5, b = 5),
    "nrow\\(X\\) must equal length\\(Y\\)"
  )
})

# ============================================================
# 15. Input validation: invalid alpha
# ============================================================
test_that("compute_bscb_hmc errors on alpha outside (0, 1)", {
  d <- make_test_data(p = 2)
  expect_error(
    compute_bscb_hmc(Y = d$Y, X = d$X, V = diag(d$n),
                     a = -5, b = 5, alpha = 1.5),
    "alpha must be in"
  )
})

# ============================================================
# 16. Input validation: a >= b
# ============================================================
test_that("compute_bscb_hmc errors when a >= b", {
  d <- make_test_data(p = 2)
  expect_error(
    compute_bscb_hmc(Y = d$Y, X = d$X, V = diag(d$n), a = 5, b = -5),
    "a must be less than b"
  )
})

# ============================================================
# 17. params list stores all settings correctly
# ============================================================
test_that("params list stores all HMC settings", {
  skip_on_cran()
  skip_if_not(instantiate::stan_cmdstan_exists(), "CmdStan not available")

  d   <- make_test_data(p = 2)
  fit <- do.call(compute_bscb_hmc, c(
    list(Y = d$Y, X = d$X, V = diag(d$n), a = -5, b = 5,
         prior_type = "normal_half_cauchy", AR_setting = 0),
    HMC_fast
  ))

  expect_equal(fit$params$prior_type,    "normal_half_cauchy")
  expect_equal(fit$params$AR_setting,    0)
  expect_equal(fit$params$L,             HMC_fast$L)
  expect_equal(fit$params$draw_num,      HMC_fast$draw_num)
  expect_equal(fit$params$iter_sampling, HMC_fast$iter_sampling)
  expect_equal(fit$params$chains,        HMC_fast$chains)
})

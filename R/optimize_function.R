##########################################################################################
####################################.     Functions   ####################################
##########################################################################################

# =============================================================================
# Optimization and helper functions for BSCB/BPCB computation
# =============================================================================

#' Create a polynomial basis vector function
#'
#' Returns a function that maps a scalar \eqn{x} to the polynomial basis
#' vector \eqn{(1, x, x^2, \ldots, x^p)}.
#'
#' @param p Non-negative integer. Polynomial degree.
#'
#' @return A function \code{f(x)} that returns \eqn{(1, x, \ldots, x^p)}
#'   as a numeric vector (for scalar \code{x}) or matrix (for vector \code{x}).
#'
#' @export
#'
#' @examples
#' f <- create_order_form(p = 2)
#' f(3)   # returns c(1, 3, 9)
#' f(c(1, 2, 3))  # returns a matrix
create_order_form <- function(p) {
  if (!is.numeric(p) || p < 0 || p != floor(p)) {
    stop("p must be a non-negative integer")
  }

  if (p == 1) {
    function(x) {
      if (is.matrix(x) || length(x) > 1) {
        return(cbind(1, x))
      } else {
        return(c(1, x))
      }
    }
  } else if (p == 2) {
    function(x) {
      if (is.matrix(x) || length(x) > 1) {
        return(cbind(1, x, x^2))
      } else {
        return(c(1, x, x^2))
      }
    }
  } else if (p == 3) {
    function(x) {
      if (is.matrix(x) || length(x) > 1) {
        return(cbind(1, x, x^2, x^3))
      } else {
        return(c(1, x, x^2, x^3))
      }
    }
  } else {
    powers <- 0:p
    function(x) {
      result <- outer(x, powers, `^`)
      if (length(x) == 1) return(as.vector(result))
      return(result)
    }
  }
}

# -----------------------------------------------------------------------------
# T(x) statistic functions
# These functions compute the standardised test statistic T(x) used in the
# critical constant estimation and coverage evaluation.
# -----------------------------------------------------------------------------

#' T(x) for computing lambda (PSCP): uses posterior draw theta_hat
#' @param x Numeric. Evaluation point.
#' @param theta_hat Numeric vector. Posterior draw of theta.
#' @param mu_star Numeric vector. Posterior mean of theta.
#' @param cov_theta Numeric matrix. Posterior covariance of theta.
#' @return Numeric scalar. Value of T(x).
#' @export

fn_Bayes_PCP <- function(x, theta_hat, mu_star, cov_theta) {
  x_i <- order_form(x)
  numerator <- abs((x_i)%*%t(theta_hat-t(mu_star)))
  denominator <- sqrt(x_i%*%cov_theta%*%t(t(x_i)))
  lambda <- numerator%*%solve(denominator)
  return(lambda)
}

#' Negative T(x) for minimisation-based optimisers (e.g. DEoptim)
#' @inheritParams fn_Bayes_PCP
#' @return Numeric scalar. Negative value of T(x).
#' @export
fn_neg_Bayes_PCP <- function(x, theta_hat, mu_star, cov_theta) {
  x_i <- order_form(x)
  numerator <- abs((x_i)%*%t(theta_hat-t(mu_star)))
  denominator <- sqrt(x_i%*%cov_theta%*%t(t(x_i)))
  lambda <- numerator%*%solve(denominator)*(-1)
  return(lambda)
}

#' T(x) for computing ESCR: uses true parameter theta_true
#' @param x Numeric. Evaluation point.
#' @param theta_true Numeric vector. True regression coefficients.
#' @param mu_star Numeric vector. Posterior mean of theta.
#' @param cov_theta Numeric matrix. Posterior covariance of theta.
#' @return Numeric scalar. Value of T(x).
#' @export
fn_Bayes_ECR <- function(x, theta_true, mu_star, cov_theta) {
  x_i <- order_form(x)
  numerator <- abs((x_i)%*%t(theta_true-t(mu_star)))
  denominator <- sqrt(x_i%*%cov_theta%*%t(t(x_i)))
  lambda <- numerator%*%solve(denominator)
  return(lambda)
}

#' T(x) for computing ESCR for DEoptim: uses true parameter theta_true
#' @param x Numeric. Evaluation point.
#' @param theta_true Numeric vector. True regression coefficients.
#' @param mu_star Numeric vector. Posterior mean of theta.
#' @param cov_theta Numeric matrix. Posterior covariance of theta.
#' @return Numeric scalar. Value of T(x).
#' @export
fn_neg_Bayes_ECR <- function(x, theta_true, mu_star, cov_theta) {
  x_i <- order_form(x)
  numerator <- abs((x_i)%*%t(theta_true-t(mu_star)))
  denominator <- sqrt(x_i%*%cov_theta%*%t(t(x_i)))
  lambda <- numerator%*%solve(denominator)*(-1)
  return(lambda)
}


#' T(x) to compute ESCR for frequentist methods
#' @param x Numeric. Evaluation point.
#' @param theta_true Numeric vector. True regression coefficients.
#' @param lm_theta_hat Numeric vector. OLS estimate of theta.
#' @param S Numeric. Residual standard error.
#' @param inv Numeric matrix. Inverse of \eqn{X^T X}.
#' @return Numeric scalar. Value of T(x).
#' @export
fn_Freq_ECR <- function(x, theta_true, lm_theta_hat, S, inv) {
  x_i <- order_form(x)
  numerator <- abs((x_i) %*% t(theta_true - t(lm_theta_hat)))
  denominator <- sqrt(x_i%*%(as.numeric(S^2) * inv)%*%t(t(x_i))) # S, inv
  lambda <- numerator %*% solve(denominator)
  return(lambda)
}

#### find the global maximum of T(x) form function
# -----------------------------------------------------------------------------
# Global maximum finders
# Two methods:
#   find_global_maximum     : grid search + local optimisation (fallback)
#   find_global_maximum_h_all: analytic polyroot method (recommended, Liu's method)
# -----------------------------------------------------------------------------

#' Find the global maximum of T(x) via grid search and local optimisation
#'
#' Fallback method using a coarse grid search combined with \code{uniroot}
#' and \code{optimize}. For most cases, \code{find_global_maximum_h_all}
#' (Liu's analytic method) is preferred.
#'
#' @param fn Function. The objective function T(x) to maximise.
#' @param a Numeric. Left endpoint of the search interval.
#' @param b Numeric. Right endpoint of the search interval.
#' @param order_form Function. Polynomial basis function from
#'   \code{create_order_form}.
#' @param theta Numeric vector. Posterior draw of theta.
#' @param mu_star Numeric vector. Posterior mean of theta.
#' @param cov_mat Numeric matrix. Covariance matrix.
#' @param tol Numeric. Numerical tolerance. Default \code{1e-6}.
#' @param n_grid Integer. Number of grid points. Default \code{100}.
#' @return A list with components \code{maximum}, \code{x_max}, and
#'   \code{all_candidates}.
#' @export
find_global_maximum <- function(fn, a, b, order_form, theta, mu_star, cov_mat,
                                tol = 1e-6, n_grid = 100) {

  numerator_fn <- function(x) {
    x_i <- order_form(x)
    as.numeric((x_i) %*% t(theta - t(mu_star)))
  }

  denominator_fn <- function(x) {
    x_i <- order_form(x)
    as.numeric(sqrt(x_i %*% cov_mat %*% t(t(x_i))))
  }

  candidate_x <- c()
  candidate_values <- c()

  candidate_x <- c(candidate_x, a, b)
  candidate_values <- c(candidate_values, fn(a), fn(b))

  grid <- seq(a, b, length.out = n_grid)
  num_values <- sapply(grid, numerator_fn)

  sign_changes <- which(diff(sign(num_values)) != 0)

  for (idx in sign_changes) {
    tryCatch({
      zero_point <- uniroot(numerator_fn,
                            interval = c(grid[idx], grid[idx + 1]),
                            tol = tol)$root
      candidate_x <- c(candidate_x, zero_point)
      candidate_values <- c(candidate_values, fn(zero_point))
    }, error = function(e) {})
  }

  fn_deriv <- function(x, h = 1e-6) {
    (fn(x + h) - fn(x - h)) / (2 * h)
  }

  search_intervals <- list()

  zero_points <- candidate_x[candidate_x > a & candidate_x < b]
  zero_points <- sort(zero_points)

  if (length(zero_points) == 0) {
    search_intervals <- list(c(a, b))
  } else {
    all_points <- c(a, zero_points, b)
    for (i in 1:(length(all_points) - 1)) {
      left <- all_points[i] + tol
      right <- all_points[i + 1] - tol
      if (right > left) {
        search_intervals[[length(search_intervals) + 1]] <- c(left, right)
      }
    }
  }

  for (interval in search_intervals) {
    tryCatch({
      result <- optimize(fn, interval = interval, maximum = TRUE)

      if (abs(fn_deriv(result$maximum)) < 0.01) {
        candidate_x <- c(candidate_x, result$maximum)
        candidate_values <- c(candidate_values, result$objective)
      }
    }, error = function(e) {})
  }

  if (length(candidate_x) > 1) {
    unique_indices <- c()
    for (i in 1:length(candidate_x)) {
      if (!any(abs(candidate_x[unique_indices] - candidate_x[i]) < tol)) {
        unique_indices <- c(unique_indices, i)
      }
    }
    candidate_x <- candidate_x[unique_indices]
    candidate_values <- candidate_values[unique_indices]
  }

  max_idx <- which.max(candidate_values)

  return(list(
    maximum = candidate_values[max_idx],
    x_max = candidate_x[max_idx],
    all_candidates = data.frame(x = candidate_x, value = candidate_values)
  ))
}


# -----------------------------------------------------------------------------
# Band evaluation functions
# f_L_SCB / f_U_SCB: vertical distance from true curve to band boundary.
# Positive value means the true curve is inside the band at that x.
# L_SCB / U_SCB: raw lower and upper band values (used in coverage integration).
# -----------------------------------------------------------------------------

#' Vertical distance from true curve to lower band boundary
#' @param x Numeric. Evaluation point.
#' @param cov_theta Numeric matrix. Posterior covariance of theta.
#' @param mu_star Numeric vector. Posterior mean of theta.
#' @param lambda_best_optim Numeric. Critical constant lambda.
#' @param theta_true Numeric vector. True regression coefficients.
#' @return Numeric. Positive if true curve is above lower bound.
#' @export
f_L_SCB <- function(x, cov_theta, mu_star, lambda_best_optim, theta_true){
  x_i <- order_form(x)
  denominator <- sqrt(x_i%*%cov_theta%*%t(t(x_i)))
  lower_bound <- x_i%*%mu_star - lambda_best_optim*denominator
  xTtheta <- x_i%*%t(t(theta_true)) # the true value
  #xTtheta <- x_i%*%t(t(theta_draw)) # the draw theta
  y_l <- xTtheta - lower_bound
  return(y_l)
}

#' Vertical distance from true curve to upper band boundary
#' @inheritParams f_L_SCB
#' @return Numeric. Positive if true curve is below upper bound.
#' @export
f_U_SCB <- function(x, cov_theta, mu_star, lambda_best_optim, theta_true){
  x_i <- order_form(x)
  denominator <- sqrt(x_i%*%cov_theta%*%t(t(x_i)))
  upper_bound <- x_i%*%mu_star + lambda_best_optim*denominator
  xTtheta <- x_i%*%t(t(theta_true)) # the true value
  #xTtheta <- x_i%*%t(t(theta_draw)) # the draw theta
  y_u <- upper_bound - xTtheta
  return(y_u)
}


#' Evaluate lower band at x
#' @param x Numeric. Evaluation point.
#' @param cov_theta Numeric matrix. Posterior covariance of theta.
#' @param mu_star Numeric vector. Posterior mean of theta.
#' @param lambda_best_optim Numeric. Critical constant lambda.
#' @return Numeric. Lower band value at x.
#' @export
L_SCB <- function(x, cov_theta, mu_star, lambda_best_optim){
  x_i <- order_form(x)
  denominator <- sqrt(x_i%*%cov_theta%*%t(t(x_i)))
  lower_bound <- x_i%*%mu_star - lambda_best_optim*denominator
  return(lower_bound)
}

#' Evaluate upper band at x
#' @inheritParams L_SCB
#' @return Numeric. Upper band value at x.
#' @export
U_SCB <- function(x, cov_theta, mu_star, lambda_best_optim){
  x_i <- order_form(x)
  denominator <- sqrt(x_i%*%cov_theta%*%t(t(x_i)))
  upper_bound <- x_i%*%mu_star + lambda_best_optim*denominator
  return(upper_bound)
}





# -----------------------------------------------------------------------------
# Analytic sup T(x) via polyroot (Liu's method, recommended)
# -----------------------------------------------------------------------------

#' Find the global maximum of T(x) analytically via polyroot (Liu's method)
#'
#' Computes \eqn{\sup_{x \in [a,b]} T(x)} by finding the stationary points of
#' \eqn{h(x) = T(x)^2 = N(x)/D(x)} via \code{polyroot}, where \eqn{N(x)} and
#' \eqn{D(x)} are polynomials. This is the recommended method.
#'
#' @param a Numeric. Left endpoint of the search interval.
#' @param b Numeric. Right endpoint of the search interval.
#' @param d Numeric vector of length \eqn{p+1}. Direction vector
#'   (\eqn{\theta - \mu^*}).
#' @param cov_mat Numeric matrix of dimension \eqn{(p+1) \times (p+1)}.
#'   Covariance matrix.
#' @return A list with components \code{maximum}, \code{x_max}, and
#'   \code{all_candidates}.
#' @export
find_global_maximum_h_all <- function(a, b, d, cov_mat) {

  # d is a (p+1) vector
  C <- cov_mat                       # (p+1)×(p+1) matrix C is the covariance matrix
  p <- length(d) - 1                 # polynomial order，quadratic: p=2, cubic: p=3

  # polynomial multiplication
  poly_mult <- function(p, q) { # p: order of the polynomial; q: order of another polynomial
    n <- length(p) + length(q) - 1
    result <- numeric(n)
    for (i in seq_along(p))
      result[i:(i + length(q) - 1)] <- result[i:(i + length(q) - 1)] + p[i] * q
    result
  }

  # polynomial derivatives
  poly_deriv <- function(coef) {
    n <- length(coef)
    if (n == 1) return(0)
    coef[-1] * seq_len(n - 1) # coef[-1]: delete the constant term
  }

  # N(x) = (d^T x)^2 = (sum d[i] x^(i-1))^2
  lin_coef <- d
  N_coef <- poly_mult(lin_coef, lin_coef)

  # D(x) = x^T C x = sum_{i,j} C[i,j] x^(i+j-2)
  D_coef <- numeric(2 * p + 1)
  for (i in 1:(p + 1))
    for (j in 1:(p + 1))
      D_coef[i + j - 1] <- D_coef[i + j - 1] + C[i, j]

  # h'(x) = 0  <=>  N'D - ND' = 0
  Np <- poly_deriv(N_coef)
  Dp <- poly_deriv(D_coef)

  NpD <- poly_mult(Np, D_coef)
  NDp <- poly_mult(N_coef, Dp)

  len <- max(length(NpD), length(NDp))
  NpD <- c(NpD, rep(0, len - length(NpD)))
  NDp <- c(NDp, rep(0, len - length(NDp)))
  poly_eq <- NpD - NDp

  # find real roots within (a,b)
  roots <- polyroot(poly_eq)
  real_roots <- Re(roots[abs(Im(roots)) < 1e-6])
  interior_roots <- real_roots[real_roots > a & real_roots < b]


  candidates <- c(a, b, interior_roots)

  # compute T(x) = sqrt(N(x)/D(x))
  T_func <- function(x) {
    xvec <- x^(0:p)
    num  <- as.numeric(d %*% xvec)^2
    den  <- as.numeric(t(xvec) %*% C %*% xvec)
    sqrt(num / den)
  }

  values  <- sapply(candidates, T_func)
  max_idx <- which.max(values)

  return(list(
    maximum        = values[max_idx],
    x_max          = candidates[max_idx],
    all_candidates = data.frame(x = candidates, value = values)
  ))
}


#' To compute the critical constant for BSCB and BPCB; To compute PSCP for BSCB and BPCB
#' @param a Numeric. Left endpoint.
#' @param b Numeric. Right endpoint.
#' @param theta_hat Numeric vector. Posterior draw of theta.
#' @param mu_star Numeric vector. Posterior mean of theta.
#' @param cov_mat Numeric matrix. Covariance matrix.
#' @return List with \code{maximum} and \code{x_max}.
#' @export
sup_T_Bayes_PSCP <- function(a, b, theta_hat, mu_star, cov_mat) {
  d <- theta_hat - as.numeric(mu_star)
  find_global_maximum_h_all(a, b, d = d, cov_mat = cov_mat)
}

#' To compute ESCR for BSCB and BPCB
#' For BSCB use \code{cov_mat = cov_theta};
#' for BPCB use \code{cov_mat = scale_mat}.
#'
#' @param a Numeric. Left endpoint.
#' @param b Numeric. Right endpoint.
#' @param theta_true Numeric vector. True regression coefficients.
#' @param mu_star Numeric vector. Posterior mean of theta.
#' @param cov_mat Numeric matrix. Covariance matrix.
#' @return List with \code{maximum} and \code{x_max}.
#' @export
sup_T_Bayes_ESCR <- function(a, b, theta_true, mu_star, cov_mat) {
  d <- theta_true - as.numeric(mu_star)
  find_global_maximum_h_all(a, b, d = d, cov_mat = cov_mat)
}



#' To compute the critical constant for simFSCB
#' @param a Numeric. Left endpoint.
#' @param b Numeric. Right endpoint.
#' @param W_sample Numeric vector. Simulated draw.
#' @param cov_mat Numeric matrix. Inverse of \eqn{X^T X}.
#' @return List with \code{maximum} and \code{x_max}.
#' @export
sup_T_simFSCB <- function(a, b, W_sample, cov_mat) { # cov_mat <- XtX_inv
  d <- as.numeric(W_sample)
  find_global_maximum_h_all(a, b, d = d, cov_mat = cov_mat)
}

#' To compute the ESCR for FSCB and FPCB
#' @param a Numeric. Left endpoint.
#' @param b Numeric. Right endpoint.
#' @param theta_true Numeric vector. True regression coefficients.
#' @param lm_theta_hat Numeric vector. OLS estimate of theta.
#' @param cov_mat Numeric matrix. Scaled covariance matrix
#'   (\eqn{S^2 \times (X^T X)^{-1}}).
#' @return List with \code{maximum} and \code{x_max}.
#' @export
sup_T_Freq_ESCR <- function(a, b, theta_true, lm_theta_hat, cov_mat){ # cov_mat <- (as.numeric(S^2) * inv
  d <- theta_true - as.numeric(lm_theta_hat)
  find_global_maximum_h_all(a, b, d = d, cov_mat = cov_mat)
}




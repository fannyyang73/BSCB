#' BSCB-C: the Bayesian simultaneous credible band using the normal-gamma conjugate prior
#'
#' @param X
#' @param Y
#' @param alpha
#' @param a
#' @param b
#' @param L
#' @param AR_setting
#' @param rho
#' @param mu_setting
#' @param P_setting
#' @param c_value
#' @param alpha_0
#' @param optimize_type
#' @param verbose
#'
#' @return
#' @export
#'
#' @examples
compute_bscb_conjugate <- function(X, # X is a n\times (p+1) matrix
                                   Y, # Y is a n dimensional vector
                                   alpha = 0.05,
                                   a = NULL,
                                   b = NULL,
                                   L = 50000,
                                   AR_setting = 0, # 0: iid error; 1: autoregressive error
                                   rho = NULL,
                                   mu_setting = 1,
                                   P_setting = 1,
                                   c_value = 0.001,
                                   alpha_0 = 1,
                                   optimize_type = "D", # D: Doptimize; G: Global-optimize
                                   verbose = TRUE
                                   ){

  # ============ Input validation ============
  if (!is.matrix(X)) X <- as.matrix(X)
  if (!is.numeric(Y)) Y <- as.numeric(Y)
  n <- nrow(X)
  p <- ncol(X) - 1
  # ============ Define order form function ============

  order_form <- create_order_form(p)

  if (length(Y) != n) {
    stop("Length of Y must equal number of rows in X")
  }

  # Infer x range from data if not provided
  if (is.null(a)) a <- min(X[, 2])
  if (is.null(b)) b <- max(X[, 2])

  # Check AR setting
  if (AR_setting == 1 && is.null(rho)) {
    stop("rho must be provided when AR_setting = 1")
  }

  # ============ Setup covariance matrix ============
  if (AR_setting == 0) {
    V <- diag(n)
  } else if (AR_setting == 1) {
    V <- rho^abs(outer(1:n, 1:n, "-")) / (1 - rho^2)
  } else {
    stop("AR_setting must be 0 or 1")
  }
  V_inv <- solve(V)

  # ============ Compute posterior parameters ============
  XtVX <- t(X) %*% V_inv %*% X
  XtVX_inv <- solve(XtVX)

  # Estimate mu based on mu_setting
  if (mu_setting == 1) {
    # OLS
    mu <- XtVX_inv %*% t(X) %*% V_inv %*% Y
  } else if (mu_setting == 2) {
    # Huber regression
    if (!requireNamespace("MASS", quietly = TRUE)) {
      stop("Package 'MASS' required for Huber regression")
    }
    k_vals <- seq(1, 3, length.out = 10)
    cv_result <- cv_huber(X[, -1, drop = FALSE], Y, k_vals)
    huber_model <- MASS::rlm(Y ~ X[, -1], psi = MASS::psi.bisquare)
    mu <- as.matrix(coef(huber_model))
  } else if (mu_setting == 3) {
    # Ridge regression
    if (!requireNamespace("glmnet", quietly = TRUE)) {
      stop("Package 'glmnet' required for Ridge regression")
    }
    cv_ridge <- glmnet::cv.glmnet(
      x = X[, -1, drop = FALSE],
      y = Y,
      alpha = 0,
      lambda = seq(0.01, 1.5, 0.01)
    )
    mu <- as.matrix(coef(cv_ridge, s = "lambda.min"))
  } else {
    stop("mu_setting must be 1, 2, or 3")
  }

  # Compute variance and precision matrix
  error <- Y - X %*% mu
  variance <- as.numeric(t(error) %*% V_inv %*% error / (n - p - 1))

  if (P_setting == 1) {
    cal_P <- c_value * diag(p + 1)
  } else if (P_setting == 2) {
    l <- variance / sum(diag(XtVX))
    cal_P <- solve(XtVX + l * diag(p + 1))
  } else {
    stop("P_setting must be 1 or 2")
  }


  # Posterior mean and covariance
  part1_1 <- solve(XtVX + cal_P)
  part1_2 <- t(X) %*% V_inv %*% Y + cal_P %*% mu
  mu_star <- part1_1 %*% part1_2

  beta_0 <- variance
  part2_1 <- (XtVX + cal_P) * (n + 2 * alpha_0)
  part2_2 <- 2 * beta_0 + t(Y) %*% V_inv %*% Y + t(mu) %*% cal_P %*% mu - t(part1_2) %*% part1_1 %*% part1_2
  the_scalar <- 1 / part2_2
  D_star <- part2_1 * drop(the_scalar)

  dof <- n + 2 * alpha_0
  the_value <- dof / (dof - 2)
  cov_theta <- drop(the_value) * solve(D_star)
  scale_mat <- solve(D_star)

  # ============ Compute lambda via Monte Carlo ============
  if (verbose) message("Computing lambda via Monte Carlo sampling...")

  lambda_samples <- numeric(L)

  fn_Bayes_PCP <- function(x, theta_hat, mu_star, cov_theta) {
    x_i <- order_form(x)
    numerator <- abs((x_i)%*%t(theta_hat-t(mu_star)))
    denominator <- sqrt(x_i%*%cov_theta%*%t(t(x_i)))
    lambda <- numerator%*%solve(denominator)
    return(lambda)
  }

  fn_neg_Bayes_PCP <- function(x, theta_hat, mu_star, cov_theta) {
    x_i <- order_form(x)
    numerator <- abs((x_i)%*%t(theta_hat-t(mu_star)))
    denominator <- sqrt(x_i%*%cov_theta%*%t(t(x_i)))
    lambda <- numerator%*%solve(denominator)*(-1)
    return(lambda)
  }


  for (j in 1:L) {
    theta_hat <- mvtnorm::rmvt(n = 1, delta = mu_star, sigma = scale_mat, df = dof, type = "shifted")
    if(optimize_type == "D"){
      result <- DEoptim::DEoptim(fn = function(x) fn_neg_Bayes_PCP(x, theta_hat, mu_star, cov_theta),
                        lower = a,
                        upper = b)
      neg_optim <- result$optim$bestval
      lambda_samples[j] <- neg_optim*(-1)

    }else if(optimize_type == "G"){
      result <- find_global_maximum(fn = function(x) fn_Bayes_PCP(x, theta_hat, mu_star, cov_theta),
                                    a, b, order_form, theta = theta_hat, mu_star = mu_star, cov_mat = cov_theta)
      lambda_samples[j] <- result$maximum
    }
  }
  lambda <- quantile(lambda_samples, probs = 1 - alpha)
  if (verbose) {
    message("The critical constant lambda = ", round(lambda, 6))
  }

  # ============ Create bound functions ============
  lower_bound <- function(x) {
    if (length(x) == 1) {
      x_i <- matrix(order_form(x), ncol = 1)
      std_error <- sqrt(t(x_i) %*% cov_theta %*% x_i)
      return(as.numeric(t(x_i) %*% mu_star - lambda * std_error))
    } else {
      # Vectorized version
      sapply(x, function(xi) {
        x_i <- matrix(order_form(xi), ncol = 1)
        std_error <- sqrt(t(x_i) %*% cov_theta %*% x_i)
        return(as.numeric(t(x_i) %*% mu_star - lambda * std_error))
      })
    }
  }

  upper_bound <- function(x) {
    if (length(x) == 1) {
      x_i <- matrix(order_form(x), ncol = 1)
      std_error <- sqrt(t(x_i) %*% cov_theta %*% x_i)
      return(as.numeric(t(x_i) %*% mu_star + lambda * std_error))
    } else {
      # Vectorized version
      sapply(x, function(xi) {
        x_i <- matrix(order_form(xi), ncol = 1)
        std_error <- sqrt(t(x_i) %*% cov_theta %*% x_i)
        return(as.numeric(t(x_i) %*% mu_star + lambda * std_error))
      })
    }
  }

  # ============ Return bscb_fit object ============
  result <- structure(
    list(
      # Main outputs
      lambda = as.numeric(lambda),
      lower_bound = lower_bound,
      upper_bound = upper_bound,

      # Posterior parameters
      mu_star = as.vector(mu_star),
      cov_theta = cov_theta,
      dof = dof,

      # Data range
      x_range = c(a, b),

      # Metadata
      call = match.call(),
      method = "conjugate",
      n = n,
      p = p,
      alpha = alpha,

      # Data (for plotting)
      data = list(X = X, Y = Y),

      # Additional info
      lambda_samples = lambda_samples,
      params = list(
        AR_setting = AR_setting,
        rho = rho,
        mu_setting = mu_setting,
        P_setting = P_setting,
        c_value = c_value,
        alpha_0 = alpha_0,
        L = L,
        optimize_type = optimize_type
      )
    ),
    class = "bscb_fit"
  )

  return(result)
}




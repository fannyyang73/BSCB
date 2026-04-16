#' BSCB-J: Bayesian Simultaneous Credible Band under the Independent Jeffreys Prior
#'
#' Constructs a \eqn{(1 - \alpha)} two-sided Bayesian simultaneous credible band
#' for polynomial regression using the independent Jeffreys prior. The marginal
#' posterior of \eqn{\theta} follows a multivariate-t distribution with degrees
#' of freedom \eqn{n - p - 1}.
#'
#' @param X Numeric matrix of dimension \eqn{n \times (p+1)}. Design matrix
#'   with intercept in the first column.
#' @param Y Numeric vector of length \eqn{n}. Response variable.
#' @param alpha Numeric. Nominal mis-coverage level; the band targets
#'   \eqn{1 - \alpha} simultaneous coverage. Default is \code{0.05}.
#' @param a Numeric. Left endpoint of the covariate domain \eqn{[a, b]}.
#'   Inferred from \code{X[, 2]} if \code{NULL}.
#' @param b Numeric. Right endpoint of the covariate domain \eqn{[a, b]}.
#'   Inferred from \code{X[, 2]} if \code{NULL}.
#' @param L Integer. Number of Monte Carlo draws for computing the critical
#'   constant \eqn{\lambda}. Default is \code{500000}.
#' @param AR_setting Integer. Error covariance structure:
#'   \code{0} = i.i.d. errors (default);
#'   \code{1} = AR(1) errors.
#' @param rho Numeric. AR(1) coefficient. Required when \code{AR_setting = 1}.
#' @param optimize_type Character. Method for computing
#'   \eqn{\sup_{x \in [a,b]} T(x)}:
#'   \code{"P"} = polyroot analytical method (recommended);
#'   \code{"G"} = global optimisation;
#'   \code{"D"} = differential evolution (DEoptim).
#' @param theta_true Numeric vector of length \eqn{p + 1}. True regression
#'   coefficients. Optional; stored in the output for diagnostic use.
#' @param verbose Logical. If \code{TRUE} (default), prints progress messages.
#'
#' @return An object of class \code{"bscb_fit"}, a list containing:
#' \describe{
#'   \item{lambda}{Critical constant for the credible band.}
#'   \item{lower_bound}{Function: computes the lower band at a given \code{x}.}
#'   \item{upper_bound}{Function: computes the upper band at a given \code{x}.}
#'   \item{mu_star}{Posterior mean of \eqn{\theta} (GLS estimate).}
#'   \item{dof}{Degrees of freedom of the marginal posterior (\eqn{n - p - 1}).}
#'   \item{scale_mat}{Scale matrix \eqn{\Sigma_0} of the marginal
#'     multivariate-t posterior distribution of \eqn{\theta}.}
#'   \item{cov_theta}{Posterior covariance matrix of \eqn{\theta}. The posterior
#'     covariance matrix equals \eqn{\text{Cov}(\theta)=\frac{\nu}{\nu-2} \Sigma_0},
#'     where \eqn{\nu} is the degrees of freedom (\code{dof}).}
#'   \item{x_range}{Covariate domain \eqn{[a, b]}.}
#'   \item{lambda_samples}{Monte Carlo samples used to compute \eqn{\lambda}.}
#'   \item{theta_true}{True parameters (if supplied).}
#'   \item{method}{Character string \code{"independent_jeffreys"}.}
#'   \item{params}{List of configuration parameters.}
#' }
#'
#' @export
#'
#' @examples
#' \donttest{
#' # Quadratic model with i.i.d. errors
#' # This is for a quick demonstration;
#' # For actual use, please set L = 500,000.
#' set.seed(123)
#' n <- 50
#' x <- seq(-5, 5, length.out = n)
#' X <- cbind(1, x, x^2)
#' theta_true <- c(-6, -3, 0.25)
#' Y <- X %*% theta_true + rnorm(n, sd = 0.2)
#'
#' fit <- compute_bscb_ind_jeffreys(
#'   X          = X,
#'   Y          = Y,
#'   alpha      = 0.05,
#'   a          = -5,
#'   b          =  5,
#'   L          = 5000,
#'   theta_true = theta_true,
#'   verbose    = FALSE
#' )
#'
#' # Critical constant
#' fit$lambda
#'
#' # Evaluate the band over a grid
#' x_seq     <- seq(-5, 5, length.out = 200)
#' lower_vec <- fit$lower_bound(x_seq)
#' upper_vec <- fit$upper_bound(x_seq)
#' }
#'
#' \donttest{
#' # Full example with recommended L
#' fit_full <- compute_bscb_ind_jeffreys(
#'   X = X, Y = Y, alpha = 0.05, a = -5, b = 5,
#'   L = 50000, theta_true = theta_true
#' )
#' }
compute_bscb_ind_jeffreys <- function(X, # X is a n\times (p+1) matrix
                                   Y, # Y is a n dimensional vector
                                   alpha = 0.05,
                                   a = NULL,
                                   b = NULL,
                                   L = 500000,
                                   AR_setting = 0, # 0: iid error; 1: autoregressive error
                                   rho = NULL,
                                   optimize_type = c("P","G","D"), # P: Polyroot; G: Global-optimize; D: Doptimize;
                                   theta_true = NULL,
                                   verbose = TRUE
){

  optimize_type <- match.arg(optimize_type)
  # ============ Input validation ============
  if (!is.matrix(X)) X <- as.matrix(X)
  if (!is.numeric(Y)) Y <- as.numeric(Y)
  n <- nrow(X)
  p <- ncol(X) - 1
  q <- ncol(X)          # q = p + 1 (intercept included)
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
  IJ_param <- compute_IJ_param(X = X, Y = Y, V = V)

  mu_star   <- IJ_param$marginal_pos_theta_mean    # marginal posterior of theta: mean
  scale_mat <- IJ_param$marginal_pos_theta_scale   # marginal posterior of theta: scale matrix
  dof       <- IJ_param$marginal_pos_theta_dof     # marginal posterior of theta: degrees of freedom
  cov_theta <- (dof / (dof - 2)) * scale_mat       # marginal posterior of theta: covariance matrix

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
                                 upper = b,
                                 control = DEoptim::DEoptim.control(itermax = 200, NP = 50, trace = FALSE)
      )

      neg_optim <- result$optim$bestval
      lambda_samples[j] <- neg_optim*(-1)

    }else if(optimize_type == "G"){
      result <- find_global_maximum(fn = function(x) fn_Bayes_PCP(x, theta_hat, mu_star, cov_theta),
                                    a, b, order_form, theta = theta_hat, mu_star = mu_star, cov_mat = cov_theta)
      lambda_samples[j] <- result$maximum
    }else if(optimize_type == "P"){
      result <- sup_T_Bayes_PSCP(a, b, theta_hat, mu_star, cov_mat = cov_theta)
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

      theta_true = theta_true,
      order_form = order_form,

      # Posterior parameters
      mu_star = as.vector(mu_star),
      dof = dof,
      scale_mat = scale_mat,
      cov_theta = cov_theta,

      # Data range
      x_range = c(a, b),

      # Metadata
      call = match.call(),
      method = "independent_jeffreys",
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
        L = L,
        optimize_type = optimize_type
      )
    ),
    class = "bscb_fit"
  )

  return(result)
}




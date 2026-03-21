#' Compute the coverage of BSCB
#' @param fit A BSCB fit object containing lambda, mu_star, cov_theta, theta_true, x_range, order_form
#' @param verbose Logical, whether to print rejection message (default: FALSE)
#'
#' @return Integer: 1 if covered, 0 if not covered
#' @export
#'
#' @examples
#' # Setup
#' set.seed(123)
#' n <- 50
#' p <- 2
#' x <- seq(-5, 5, length.out = n)
#' X <- cbind(1, x, x^2)
#' theta_true <- c(-6, -3, 0.25)
#'
#' # Generate data and compute BSCB
#' Y <- X %*% theta_true + rnorm(n, 0, 0.2)
#' fit <- compute_bscb_conjugate(X, Y, alpha = 0.05, a = -5, b = 5,
#'                                L = 10000, theta_true = theta_true,
#'                                verbose = FALSE)
#'
#' # Check the empirical simultaneous coverage rate (ESCR)
#' is_covered <- coverage_ESCR(fit,  optimize_type ="P", verbose = TRUE)
#' cat("Coverage indicator:", is_covered, "\n")
#'
coverage_ESCR <- function(fit,
                          optimize_type = c("P","G","D"), # P: Polyroot; G: Global-optimize; D: Doptimize;
                          verbose = FALSE){

  fn_Bayes_ECR <- function(x, theta_true, mu_star, cov_theta) {
    x_i <- order_form(x)
    numerator <- abs((x_i)%*%t(theta_true-t(mu_star)))
    denominator <- sqrt(x_i%*%cov_theta%*%t(t(x_i)))
    lambda <- numerator%*%solve(denominator)
    return(lambda)
  }

  fn_neg_Bayes_ECR <- function(x, theta_true, mu_star, cov_theta) {
    x_i <- order_form(x)
    numerator <- abs((x_i)%*%t(theta_true-t(mu_star)))
    denominator <- sqrt(x_i%*%cov_theta%*%t(t(x_i)))
    lambda <- numerator%*%solve(denominator)*(-1)
    return(lambda)
  }


  lambda_best <- fit$lambda
  mu_star <- fit$mu_star
  cov_theta <- fit$cov_theta
  theta_true <- fit$theta_true
  a <- fit$x_range[1]
  b <- fit$x_range[2]
  order_form <- fit$order_form


  result_ECR <- find_global_maximum(
    fn = function(x) fn_Bayes_ECR(x, theta_true, mu_star, cov_theta),
    a = a,
    b = b,
    order_form = order_form,
    theta = theta_true,
    mu_star = mu_star,
    cov_mat = cov_theta
  )


  if(optimize_type == "D"){
    result <- DEoptim::DEoptim(fn = function(x) fn_neg_Bayes_ECR(x, theta_true, mu_star, cov_theta),
                               lower = a,
                               upper = b,
                               control = DEoptim::DEoptim.control(itermax = 200, NP = 50, trace = FALSE)
    )

    neg_optim <- result$optim$bestval
    result_ECR <- neg_optim*(-1)

  }else if(optimize_type == "G"){
    result <- find_global_maximum(fn = function(x) fn_Bayes_ECR(x, theta_true, mu_star, cov_theta),
                                  a, b, order_form, theta = theta_true, mu_star = mu_star, cov_mat = cov_theta)
    result_ECR <- result$maximum
  }else if(optimize_type == "P"){
    result <- sup_T_Bayes_ESCR(a, b, theta_true, mu_star, cov_mat = cov_theta)
    result_ECR <- result$maximum
  }

  # Check coverage
  coverage_flag <- ifelse(lambda_best >= result_ECR, 1, 0)

  return(coverage_flag)
}

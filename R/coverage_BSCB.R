#' Compute the coverage of BSCB
#' @param fit A BSCB fit object containing lambda, mu_star, cov_theta, theta_true, x_range, order_form
#' @param verbose Logical, whether to print rejection message (default: FALSE)
#' @param optimize_type Character. Method for computing
#'   \eqn{\sup_{x \in [a,b]} T(x)}:
#'   \code{"P"} = polyroot function (recommended);
#'   \code{"G"} = global optimisation;
#'   \code{"D"} = Doptimize function from package DEoptim.
#' @param verbose Logical. If \code{TRUE} (default), prints the value of the
#'   critical constant lambda.
#'
#' @return Integer: 1 if covered, 0 if not covered
#' @export
#'
#' @examples
#' \donttest{
#' # This is for a quick demonstration;
#' # For actual use, please set L = 500000.
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
#'                                L = 50000, theta_true = theta_true,
#'                                verbose = FALSE)
#'
#' # Check the empirical simultaneous coverage rate (ESCR)
#' is_covered <- coverage_ESCR(fit,  optimize_type ="P", verbose = TRUE)
#' cat("Coverage indicator:", is_covered, "\n")
#' }

coverage_ESCR <- function(fit,
                          optimize_type = c("P","G","D"), # P: Polyroot; G: Global-optimize; D: Doptimize;
                          verbose = FALSE){

  optimize_type <- match.arg(optimize_type)


  lambda_best <- fit$lambda
  mu_star <- fit$mu_star
  cov_theta <- fit$cov_theta
  theta_true <- fit$theta_true
  a <- fit$x_range[1]
  b <- fit$x_range[2]
  order_form <- fit$order_form
  if (is.null(theta_true)) {
    stop("fit$theta_true is NULL. Supply theta_true when fitting the model.")
  }

  # ============ Compute sup T(x) at theta_true ============

  if(optimize_type == "D"){ #fn_neg_Bayes_ECR
    result <- DEoptim::DEoptim(fn = function(x){
      x_i         <- order_form(x)
      numerator   <- abs(x_i %*% t(theta_true - t(mu_star)))
      denominator <- sqrt(x_i %*% cov_theta %*% t(t(x_i)))
      -as.numeric(numerator %*% solve(denominator))
    },
    lower = a,
    upper = b,
    control = DEoptim::DEoptim.control(itermax = 200, NP = 50, trace = FALSE)
    )

    neg_optim <- result$optim$bestval
    result_ECR <- neg_optim*(-1)

  }else if(optimize_type == "G"){ #fn_Bayes_ECR
    result <- find_global_maximum(fn = function(x){
      x_i         <- order_form(x)
      numerator   <- abs(x_i %*% t(theta_true - t(mu_star)))
      denominator <- sqrt(x_i %*% cov_theta %*% t(t(x_i)))
      as.numeric(numerator %*% solve(denominator))
    },
    a, b, order_form, theta = theta_true, mu_star = mu_star, cov_mat = cov_theta)
    result_ECR <- result$maximum
  }else if(optimize_type == "P"){
    result <- sup_T_Bayes_ESCR(a, b, theta_true, mu_star, cov_mat = cov_theta)
    result_ECR <- result$maximum
  }

  # Check coverage
  coverage_flag <- as.integer(lambda_best >= result_ECR)

  if (verbose && coverage_flag == 0) {
    message("Coverage failed: sup T(x) = ", round(result_ECR, 6),
            " > lambda = ", round(lambda_best, 6))
  }

  return(coverage_flag)
}



#' Compute the Posterior Simultaneous Coverage Probability (PSCP)
#'
#' Estimates the posterior simultaneous coverage probability (PSCP) of a
#' constructed BSCB by Monte Carlo integration over the posterior distribution
#' of \eqn{\theta}. For each posterior draw \eqn{\hat{\theta}}, the supremum
#' \eqn{\sup_{x \in [a,b]} T(x)} is computed and compared against the critical
#' constant \eqn{\lambda}. The PSCP is the proportion of draws for which
#' \eqn{\sup T(x) \leq \lambda}.
#'
#' @param fit An object of class \code{"bscb_fit"} returned by
#'   \code{\link{compute_bscb_conjugate}} or
#'   \code{\link{compute_bscb_ind_jeffreys}}.
#'   Must contain \code{lambda},
#'   \code{mu_star}, \code{cov_theta}, \code{dof}, \code{x_range}, and
#'   \code{order_form}.
#' @param draw_num Integer. Number of Monte Carlo draws for estimating PSCP.
#'   Default is \code{10000}.
#' @param optimize_type Character. Method for computing
#'   \eqn{\sup_{x \in [a,b]} T(x)}:
#'   \code{"P"} = polyroot analytical method (recommended);
#'   \code{"G"} = global optimisation;
#'   \code{"D"} = differential evolution (DEoptim).
#' @param verbose Logical. If \code{TRUE}, prints the estimated PSCP value.
#'   Default is \code{FALSE}.
#'
#' @return Numeric. Estimated posterior simultaneous coverage probability,
#'   a value in \eqn{[0, 1]}.
#'
#' @seealso \code{\link{coverage_ESCR}},
#'   \code{\link{compute_bscb_conjugate}},
#'   \code{\link{compute_bscb_ind_jeffreys}}
#'
#' @export
#'
#' @examples
#' # This is for a quick demonstration;
#' # For actual use, please set L = 500000 and draw_num = 10000.
#' set.seed(123)
#' n <- 50
#' x <- seq(-5, 5, length.out = n)
#' X <- cbind(1, x, x^2)
#' theta_true <- c(-6, -3, 0.25)
#' Y <- X %*% theta_true + rnorm(n, sd = 0.2)
#'
#' fit <- compute_bscb_conjugate(
#'   X          = X,
#'   Y          = Y,
#'   alpha      = 0.05,
#'   a          = -5,
#'   b          =  5,
#'   L          = 1000,
#'   theta_true = theta_true,
#'   verbose    = FALSE
#' )
#'
#' coverage_PSCP(fit, draw_num = 500, optimize_type = "P", verbose = TRUE)
#'
#' \donttest{
#' # Full example with recommended draw_num
#' coverage_PSCP(fit, draw_num = 10000, optimize_type = "P")
#' }
coverage_PSCP <- function(fit,
                          draw_num      = 10000,
                          optimize_type = c("P","G","D"), # P: Polyroot; G: Global-optimize; D: Doptimize;
                          verbose = FALSE){

  optimize_type <- match.arg(optimize_type)

  # ============ Extract components from fit object ============
  lambda_best <- fit$lambda
  mu_star <- fit$mu_star
  cov_theta <- fit$cov_theta
  is_HMC <- !is.null(fit$method) && fit$method == "HMC"
  if (is_HMC) {
    theta_mat <- fit$theta_mat
    n_samples <- nrow(theta_mat)
  } else {
    dof       <- fit$dof
    scale_mat <- fit$scale_mat
  }
  theta_true <- fit$theta_true
  a <- fit$x_range[1]
  b <- fit$x_range[2]
  order_form <- fit$order_form
  if (is.null(theta_true)) {
    stop("fit$theta_true is NULL. Supply theta_true when fitting the model.")
  }

  # ============ Monte Carlo estimation of PSCP ============
  # fn_Bayes_PCP and fn_neg_Bayes_PCP are defined in optimize_function.R
  cover_num <- 0
  for (j in 1:draw_num){
    if (is_HMC){
      theta_hat <- theta_mat[sample(n_samples, 1), ]
    }else{
      theta_hat <- mvtnorm::rmvt(n = 1, delta = mu_star, sigma = scale_mat, df = dof, type = "shifted")
    }


    if(optimize_type == "D"){#fn_neg_Bayes_PCP
      result <- DEoptim::DEoptim(fn = function(x){
        x_i         <- order_form(x)
        numerator   <- abs(x_i %*% t(theta_hat - t(mu_star)))
        denominator <- sqrt(x_i %*% cov_theta %*% t(t(x_i)))
        -as.numeric(numerator %*% solve(denominator))
      },
      lower = a,
      upper = b,
      control = DEoptim::DEoptim.control(itermax = 200, NP = 50, trace = FALSE)
      )

      neg_optim <- result$optim$bestval
      result_PSCP <- neg_optim*(-1)

    }else if(optimize_type == "G"){ #fn_Bayes_PCP
      result <- find_global_maximum(fn = function(x){
        x_i         <- order_form(x)
        numerator   <- abs(x_i %*% t(theta_hat - t(mu_star)))
        denominator <- sqrt(x_i %*% cov_theta %*% t(t(x_i)))
        as.numeric(numerator %*% solve(denominator))
      },
      a, b, order_form, theta = theta_hat, mu_star = mu_star, cov_mat = cov_theta)
      result_PSCP <- result$maximum
    }else if(optimize_type == "P"){
      result <- sup_T_Bayes_PSCP(a, b, theta_hat, mu_star, cov_mat = cov_theta)
      result_PSCP <- result$maximum
    }

    if (result_PSCP <= lambda_best) cover_num <- cover_num + 1
  }
  coverage_draw <- cover_num / draw_num
  if (verbose) message("PSCP = ", round(coverage_draw, 6))

  return(coverage_draw)
}

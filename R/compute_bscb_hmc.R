#' Compute BSCB via Hamiltonian Monte Carlo
#'
#' @description
#' Constructs a Bayesian Simultaneous Credible Band (BSCB) for polynomial
#' regression under non-conjugate priors using Hamiltonian Monte Carlo (HMC) via Stan.
#' Supports two prior specifications: Normal-Normal and Normal-half-Cauchy.
#' The critical constant lambda is estimated by Monte Carlo, and the
#' Posterior Simultaneous Coverage Probability (PSCP) is also returned.
#' @importFrom stats cov
#' @importFrom posterior as_draws_df
#'
#' @param Y Numeric vector of responses of length \code{n}.
#' @param X Design matrix of dimension \code{n x (p+1)}, including an
#'   intercept column. The second column must contain the raw covariate values.
#' @param V Error covariance matrix of dimension \code{n x n}.
#'   Use \code{diag(n)} for i.i.d. errors (default).
#' @param a Left endpoint of the covariate domain \code{[a, b]}.
#' @param b Right endpoint of the covariate domain \code{[a, b]}.
#' @param theta_true Numeric vector of true regression coefficients of length
#'   \code{p+1}. Used to evaluate ESCR in simulation studies. Set to
#'   \code{NULL} (default) when the true coefficients are unknown.
#' @param alpha Nominal miscoverage rate. The credible band targets
#'   \code{1 - alpha} simultaneous coverage. Default is \code{0.05}.
#' @param prior_type Character string specifying the prior on
#'   \code{(theta, sigma)}. Either \code{"normal_half_cauchy"} (default,
#'   recommended) or \code{"normal_normal"}.
#' @param normal_theta_sd Prior standard deviation for each component of
#'   \code{theta} under the Normal-Normal prior. Default is \code{10}.
#' @param normal_sigma_sd Prior standard deviation for \code{sigma} under
#'   the Normal-Normal prior. Default is \code{5}.
#' @param cauchy_scale Scale parameter of the half-Cauchy prior on
#'   \code{sigma} under the Normal-half-Cauchy prior. Default is \code{2}.
#' @param iter_sampling Number of post-warmup HMC draws per chain.
#'   Default is \code{4000}.
#' @param iter_warmup Number of warmup draws per chain. Default is \code{4000}.
#' @param chains Number of Markov chains. Default is \code{4}.
#' @param thin_number Positive integer. Thinning interval for posterior draws.
#'   A value of \code{k} retains every \code{k}-th draw from each chain.
#'   Default is \code{1} (no thinning).
#' @param adapt_delta Target acceptance probability for the NUTS sampler.
#'   Default is \code{0.95}.
#' @param max_treedepth Maximum tree depth for the NUTS sampler.
#'   Default is \code{15}.
#' @param AR_setting Integer. \code{0} for i.i.d. errors (default),
#'   \code{1} for AR(1) errors.
#' @param rho AR(1) autocorrelation coefficient. Only used when
#'   \code{AR_setting == 1}. Default is \code{0}.
#' @param optimize_type Character. Method for computing
#'   \eqn{\sup_{x \in [a,b]} T(x)}:
#'   \code{"P"} = polyroot analytical method (recommended, default);
#'   \code{"G"} = global optimisation (grid-based);
#'   \code{"D"} = differential evolution (\code{DEoptim}).
#' @param L Number of Monte Carlo draws used to estimate the critical
#'   constant \code{lambda}. Default is \code{500000}.
#' @param draw_num Number of Monte Carlo draws used to estimate the PSCP.
#'   Default is \code{10000}.
#'
#' @return An object of class \code{"bscb_fit"}, which is a list with the
#'   following components:
#'   \item{lambda}{Estimated critical constant at level \code{1 - alpha}.}
#'   \item{lower_bound}{Lower credible band evaluated on a fine grid over
#'     \code{[a, b]}.}
#'   \item{upper_bound}{Upper credible band evaluated on a fine grid over
#'     \code{[a, b]}.}
#'   \item{theta_true}{True regression coefficients (if supplied).}
#'   \item{order_form}{Polynomial order form used internally.}
#'   \item{mu_star}{Posterior mean of \code{theta} (length \code{p+1}).}
#'   \item{cov_theta}{Posterior covariance matrix of \code{theta}.}
#'   \item{theta_mat}{Matrix of posterior draws,
#'     \code{(chains * floor(iter_sampling / thin_number)) x (p+1)}.}
#'   \item{x_range}{Numeric vector \code{c(a, b)}.}
#'   \item{call}{The matched call.}
#'   \item{method}{Character string \code{"HMC"}.}
#'   \item{n}{Sample size.}
#'   \item{p}{Polynomial degree.}
#'   \item{alpha}{Nominal miscoverage rate.}
#'   \item{data}{List containing the design matrix \code{X} and response
#'     vector \code{Y}.}
#'   \item{lambda_samples}{Numeric vector of length \code{L} containing the
#'     Monte Carlo supremum draws used to derive \code{lambda}.}
#'   \item{params}{List of additional settings: \code{AR_setting},
#'     \code{rho}, \code{prior_type}, \code{normal_theta_sd},
#'     \code{normal_sigma_sd}, \code{cauchy_scale}, \code{iter_sampling},
#'     \code{iter_warmup}, \code{chains}, \code{thin_number}, \code{L},
#'     \code{draw_num}, \code{optimize_type}.}
#'
#'
#' @seealso \code{\link{compute_bscb_conjugate}}, \code{\link{compute_bscb_ind_jeffreys}}
#'
#' @examples
#' \donttest{
#'   set.seed(42)
#'   n <- 20; p <- 2
#'   x_seq <- seq(-5, 5, length.out = n)
#'   X <- cbind(1, x_seq, x_seq^2)
#'   theta_true <- c(-6, -3, 0.25)
#'   Y <- as.numeric(X %*% theta_true + rnorm(n, sd = 0.2))
#'   fit <- compute_bscb_hmc(
#'     Y = Y, X = X, V = diag(n),
#'     a = -5, b = 5,
#'     theta_true = theta_true,
#'     prior_type = "normal_half_cauchy",
#'     L = 1000, draw_num = 500   # small values for illustration only
#'   )
#'   fit$lambda
#'   fit$params$prior_type
#' }
#'
#' @export
compute_bscb_hmc <- function(Y, X, V = diag(nrow(X)),
                             a, b,
                             theta_true      = NULL,
                             alpha           = 0.05,
                             prior_type      = c("normal_half_cauchy", "normal_normal"),
                             normal_theta_sd = 10,
                             normal_sigma_sd =  5,
                             cauchy_scale    =  2,
                             iter_sampling   = 4000,
                             iter_warmup     = 4000,
                             chains          = 4,
                             thin_number     = 1,
                             adapt_delta     = 0.95,
                             max_treedepth   = 15,
                             AR_setting      = 0,
                             rho             = 0,
                             optimize_type   = c("P", "G", "D"),
                             L               = 500000,
                             draw_num        = 10000) {

  # ------------------------------------------------------------------
  # Input validation
  # ------------------------------------------------------------------
  prior_type    <- match.arg(prior_type)
  optimize_type <- match.arg(optimize_type)

  stopifnot(
    "Y must be a numeric vector"          = is.numeric(Y) && is.vector(Y),
    "X must be a numeric matrix"          = is.numeric(X) && is.matrix(X),
    "nrow(X) must equal length(Y)"        = nrow(X) == length(Y),
    "V must be a square matrix"           = is.matrix(V) && nrow(V) == ncol(V),
    "nrow(V) must equal length(Y)"        = nrow(V) == length(Y),
    "a must be less than b"               = a < b,
    "alpha must be in (0, 1)"             = alpha > 0 && alpha < 1,
    "AR_setting must be 0 or 1"           = AR_setting %in% c(0L, 1L),
    "L must be a positive integer"        = L > 0,
    "draw_num must be a positive integer" = draw_num > 0
  )

  n <- nrow(X)
  p <- ncol(X) - 1L

  stopifnot("p must be 2 or 3" = p %in% c(2L, 3L))

  if (!is.null(theta_true))
    stopifnot("theta_true must have length p+1" = length(theta_true) == p + 1L)
  if (!requireNamespace("cmdstanr", quietly = TRUE)) {
    stop(
      'Package "cmdstanr" is required for HMC methods.\n',
      'Install with: install.packages("cmdstanr", ',
      'repos = c("https://stan-dev.r-universe.dev", getOption("repos")))',
      call. = FALSE
    )
  }
  HMC_prior_type <- if (prior_type == "normal_normal") 0L else 1L

  # Store the call for the returned object
  mc <- match.call()

  # ------------------------------------------------------------------
  # Polynomial order form (used by sup_T functions and band closures)
  # ------------------------------------------------------------------
  order_form <- create_order_form(p)

  # ------------------------------------------------------------------
  # Load pre-compiled Stan model via instantiate.
  # The file bin/stan/HMC_model.stan is compiled once at install time.
  # ------------------------------------------------------------------
  mod <- withCallingHandlers(
    instantiate::stan_package_model(
      name    = "HMC_model",
      package = "BSCB"
    ),
    stan_deprecate = function(w) {
      invokeRestart("muffleWarning")
    }
  )

  # ------------------------------------------------------------------
  # Empirical Bayes hyperparameters for theta prior mean.
  # (used only when prior_type == "normal_half_cauchy")
  # ------------------------------------------------------------------
  E_NG <- compute_NG_param(
    X              = X,
    Y              = Y,
    V              = V,
    hyperparameter = "empirical"
  )

  # ------------------------------------------------------------------
  # Stan data list.
  # theta_prior_mean is a vector of length p+1, matching the unified
  # Stan model which declares  vector[p+1] theta_prior_mean.
  # ------------------------------------------------------------------
  stan_data <- list(
    p                = p,
    N                = n,
    x                = X[, 2],
    y                = as.numeric(Y),
    V                = V,
    theta_prior_sd   = sqrt(E_NG$beta_0),
    theta_prior_mean = as.numeric(E_NG$beta_hat),  # length p+1
    use_likelihood   = 1L,
    HMC_prior_type   = HMC_prior_type,
    normal_theta_sd  = normal_theta_sd,
    normal_sigma_sd  = normal_sigma_sd,
    cauchy_scale     = cauchy_scale
  )

  # ------------------------------------------------------------------
  # HMC sampling.
  # Temporary output directory is cleaned up on function exit.
  # ------------------------------------------------------------------
  output_dir <- file.path(tempdir(), paste0("bscb_hmc_", Sys.getpid()))
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(output_dir, recursive = TRUE), add = TRUE)

  fit <- mod$sample(
    data            = stan_data,
    iter_sampling   = iter_sampling,
    iter_warmup     = iter_warmup,
    chains          = chains,
    thin            = thin_number,
    parallel_chains = 1L,
    refresh         = 0L,
    adapt_delta     = adapt_delta,
    max_treedepth   = max_treedepth,
    output_dir      = output_dir
  )

  # ------------------------------------------------------------------
  # Extract posterior draws -> theta_mat: (n_samples x (p+1))
  # ------------------------------------------------------------------
  posterior_samples <- fit$draws(format = "df")
  theta_mat <- do.call(cbind, lapply(seq_len(p + 1L), function(k) {
    posterior_samples[[paste0("theta[", k, "]")]]
  }))
  n_samples <- nrow(theta_mat)
  mu_star   <- colMeans(theta_mat)
  cov_theta <- cov(theta_mat)

  # ------------------------------------------------------------------
  # Step 2: Monte Carlo estimation of critical constant lambda
  # ------------------------------------------------------------------
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

  lambda_samples <- numeric(L)
  for (j in seq_len(L)) {
    theta_hat <- theta_mat[sample(n_samples, 1L), ]
    if (optimize_type == "D") {
      res <- DEoptim::DEoptim(
        fn      = function(x) fn_neg_Bayes_PCP(x, theta_hat, mu_star, cov_theta),
        lower   = a,
        upper   = b,
        control = DEoptim::DEoptim.control(itermax = 200, NP = 50, trace = FALSE)
      )
      lambda_samples[j] <- res$optim$bestval * (-1)
    } else if (optimize_type == "G") {
      res <- find_global_maximum(
        fn = function(x) fn_Bayes_PCP(x, theta_hat, mu_star, cov_theta),
        a, b, order_form,
        theta = theta_hat, mu_star = mu_star, cov_mat = cov_theta
      )
      lambda_samples[j] <- res$maximum
    } else {  # "P"
      res <- sup_T_Bayes_PSCP(a, b, theta_hat, mu_star, cov_mat = cov_theta)
      lambda_samples[j] <- res$maximum
    }
  }
  lambda <- unname(quantile(lambda_samples, probs = 1 - alpha))

  # ------------------------------------------------------------------
  # Step 3: PSCP estimation
  # ------------------------------------------------------------------
  cover_num <- 0L
  for (j in seq_len(draw_num)) {
    theta_hat <- theta_mat[sample(n_samples, 1L), ]
    sup_val   <- if (optimize_type == "D") {
      res <- DEoptim::DEoptim(
        fn      = function(x) fn_neg_Bayes_PCP(x, theta_hat, mu_star, cov_theta),
        lower   = a,
        upper   = b,
        control = DEoptim::DEoptim.control(itermax = 200, NP = 50, trace = FALSE)
      )
      res$optim$bestval * (-1)
    } else if (optimize_type == "G") {
      find_global_maximum(
        fn = function(x) fn_Bayes_PCP(x, theta_hat, mu_star, cov_theta),
        a, b, order_form,
        theta = theta_hat, mu_star = mu_star, cov_mat = cov_theta
      )$maximum
    } else {  # "P"
      sup_T_Bayes_PSCP(a, b, theta_hat, mu_star, cov_mat = cov_theta)$maximum
    }
    if (sup_val <= lambda) cover_num <- cover_num + 1L
  }
  PSCP <- cover_num / draw_num

  # ------------------------------------------------------------------
  # Step 4: Credible band functions (closures over mu_star, cov_theta, lambda)
  # ------------------------------------------------------------------
  lower_bound <- function(x) {
    if (length(x) == 1) {
      x_i       <- matrix(order_form(x), ncol = 1)
      std_error <- sqrt(t(x_i) %*% cov_theta %*% x_i)
      return(as.numeric(t(x_i) %*% mu_star - lambda * std_error))
    } else {
      sapply(x, function(xi) {
        x_i       <- matrix(order_form(xi), ncol = 1)
        std_error <- sqrt(t(x_i) %*% cov_theta %*% x_i)
        as.numeric(t(x_i) %*% mu_star - lambda * std_error)
      })
    }
  }

  upper_bound <- function(x) {
    if (length(x) == 1) {
      x_i       <- matrix(order_form(x), ncol = 1)
      std_error <- sqrt(t(x_i) %*% cov_theta %*% x_i)
      return(as.numeric(t(x_i) %*% mu_star + lambda * std_error))
    } else {
      sapply(x, function(xi) {
        x_i       <- matrix(order_form(xi), ncol = 1)
        std_error <- sqrt(t(x_i) %*% cov_theta %*% x_i)
        as.numeric(t(x_i) %*% mu_star + lambda * std_error)
      })
    }
  }

  # ------------------------------------------------------------------
  # Return bscb_fit object
  # ------------------------------------------------------------------
  result <- structure(
    list(
      # Main outputs
      lambda      = as.numeric(lambda),
      lower_bound = lower_bound,
      upper_bound = upper_bound,
      theta_true  = theta_true,
      order_form  = order_form,
      PSCP        = PSCP,
      # Posterior parameters
      mu_star   = as.vector(mu_star),
      cov_theta = cov_theta,
      theta_mat = theta_mat,
      # Data range
      x_range = c(a, b),
      # Metadata
      call   = mc,
      method = "HMC",
      n      = n,
      p      = p,
      alpha  = alpha,
      # Data (for plotting)
      data = list(X = X, Y = Y),
      # Monte Carlo draws for lambda
      lambda_samples = lambda_samples,
      # Additional settings
      params = list(
        AR_setting      = AR_setting,
        rho             = rho,
        prior_type      = prior_type,      # "normal_half_cauchy" / "normal_normal"
        normal_theta_sd = normal_theta_sd,
        normal_sigma_sd = normal_sigma_sd,
        cauchy_scale    = cauchy_scale,
        iter_sampling   = iter_sampling,
        iter_warmup     = iter_warmup,
        chains          = chains,
        thin_number     = thin_number,
        L               = L,
        draw_num        = draw_num,
        optimize_type   = optimize_type
      )
    ),
    class = "bscb_fit"
  )

  result
}

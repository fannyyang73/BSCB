#' Generate simulation datasets for polynomial regression
#'
#' Generates a fixed design matrix X and a list of response vectors Y for
#' use in simulation studies of Bayesian simultaneous credible bands.
#' The design can be either equally-spaced (ES) or D-optimal (DO).
#'
#' @importFrom OptimalDesign od_KL
#'
#' @param p Integer. Polynomial degree. Must be 1, 2, or 3.
#' @param n Integer. Sample size.
#' @param e_sd Numeric. Error standard deviation (sigma in the paper).
#' @param theta_true Numeric vector of length \code{p + 1}. True regression
#'   coefficients.
#' @param a Numeric. Left endpoint of the covariate domain \eqn{[a, b]}.
#' @param b Numeric. Right endpoint of the covariate domain \eqn{[a, b]}.
#' @param replication Integer. Number of simulation replications.
#' @param design_index Integer. Design type:
#'   \code{1} = equally-spaced (ES);
#'   \code{2} = D-optimal (DO).
#' @param center_index Integer. Centering of covariates:
#'   \code{1} = mean-centred (default);
#'   \code{0} = uncentred;
#'   \code{2} = standardised.
#' @param n_ES_x Integer. Number of equally-spaced design points.
#'   Only used when \code{design_index = 1}.
#' @param n_DO_init_x Integer. Candidate pool size for D-optimal search.
#'   A large value (e.g. 300000) ensures that 6 support points are selected.
#'   Only used when \code{design_index = 2}.
#' @param AR_index Integer. Error structure:
#'   \code{0} = i.i.d. (default); \code{1} = AR(1).
#' @param rho Numeric. AR(1) coefficient. Only used when \code{AR_index = 1}.
#' @param batch_index Integer. Batch index used as part of the random seed
#'   (\code{set.seed(1000 * batch_index + i)} for replication \code{i}).
#' @param seed Integer or \code{NULL}. Base random seed used internally for
#' the D-optimal design search (\code{OptimalDesign::od_KL}). If \code{NULL}
#' (default), no seed is set and results may vary between runs. When
#' \code{seed} is supplied, results are fully reproducible because the
#' internal search is stopped after a fixed number of restarts
#' (\code{rest.max}), rather than after a fixed time budget.
#'
#' @return A list containing:
#' \describe{
#'   \item{X}{Design matrix of dimension \eqn{n \times (p+1)}.}
#'   \item{Y.list}{List of \code{replication} response vectors, each of
#'     length \code{n}.}
#'   \item{optimal_x}{Vector of selected support points.}
#'   \item{optimal_weights}{Vector of observation counts at each support
#'     point.}
#' }
#'
#' @export
#'
#' @examples
#'
#' \donttest{
#' # Example 1: quadratic model, D-optimal design
#' sim_data <- generate_simulation_data(
#'   p           = 2,
#'   n           = 20,
#'   e_sd        = 0.2,
#'   theta_true  = c(-6, -3, 0.25),
#'   a           = -5,
#'   b           =  5,
#'   replication = 1,
#'   design_index = 2,
#'   center_index = 1
#' )
#'
#' X      <- sim_data$X
#' Y.list <- sim_data$Y.list
#' }
#'
#' \donttest{
#' # Example 2: cubic model, equally-spaced design
#' sim_data2 <- generate_simulation_data(
#'   p            = 3,
#'   n            = 20,
#'   e_sd         = 0.2,
#'   theta_true   = c(1, 2, -1, 0.5),
#'   a            = -5,
#'   b            =  5,
#'   replication  = 1,
#'   design_index = 1,
#'   center_index = 1
#' )
#' }
generate_simulation_data <- function(p,
                                     n,
                                     e_sd,
                                     theta_true,
                                     a           = -5,
                                     b           =  5,
                                     replication = 2,
                                     design_index = 2,
                                     center_index = 1,
                                     n_ES_x       = n,
                                     n_DO_init_x  = 300000,
                                     AR_index     = 0,
                                     rho          = 0.1,
                                     batch_index  = 1,
                                     seed = NULL) {

  # ---------------------------------------------------------------------------
  # Input validation
  # ---------------------------------------------------------------------------
  stopifnot(
    "p must be 1, 2, or 3"              = p %in% 1:3,
    "design_index must be 1 or 2"       = design_index %in% 1:2,
    "center_index must be 0, 1, or 2"   = center_index %in% 0:2,
    "AR_index must be 0 or 1"           = AR_index %in% 0:1,
    "length of theta_true must be p+1"  = length(theta_true) == p + 1
  )

  # ---------------------------------------------------------------------------
  # Order form and error covariance
  # ---------------------------------------------------------------------------
  order_form <- create_order_form(p)

  if (AR_index == 0) {
    V <- diag(n)
  } else {
    V <- rho^abs(outer(1:n, 1:n, "-")) / (1 - rho^2)
  }

  # ---------------------------------------------------------------------------
  # Generate design matrix
  # ---------------------------------------------------------------------------
  if (design_index == 1) {
    # Equally-spaced design
    x_candidates <- seq(a, b, length.out = n_ES_x)
    if (center_index == 1) {
      x_candidates <- scale(x_candidates, center = TRUE, scale = FALSE)
    } else if (center_index == 2) {
      x_candidates <- scale(x_candidates, center = TRUE, scale = TRUE)
    }
    optimal_x       <- as.numeric(x_candidates)
    optimal_weights <- rep(n / n_ES_x, n_ES_x)

  } else if (design_index == 2) {
    # D-optimal design
    x_candidates <- seq(a, b, length.out = n_DO_init_x)
    if (center_index == 1) {
      x_candidates <- scale(x_candidates, center = TRUE, scale = FALSE)
    } else if (center_index == 2) {
      x_candidates <- scale(x_candidates, center = TRUE, scale = TRUE)
    }
    x_candidates <- as.numeric(x_candidates)

    design_matrix <- t(sapply(x_candidates, order_form))  # (n_DO_init_x) x (p+1)

    if(!is.null(seed)){set.seed(seed)}
    d_opt <- OptimalDesign::od_KL(
      Fx   = design_matrix,
      N    = n,
      t.max = 30,
      rest.max = 5,
      K    = 7,
      L    = 19,
      crit = "D"
    )
    optimal_x       <- x_candidates[d_opt$supp]
    optimal_weights <- d_opt$w.supp
  }

  # Build X from design
  x_rep <- rep(optimal_x, optimal_weights)
  X     <- t(sapply(x_rep, order_form))

  # ---------------------------------------------------------------------------
  # Generate response vectors
  # ---------------------------------------------------------------------------
  Y.list <- vector("list", replication)
  for (i in seq_len(replication)) {
    if(!is.null(batch_index)){set.seed(1000 * batch_index + i)}
    epsilon   <- MASS::mvrnorm(n = 1, mu = rep(0, n), Sigma = e_sd^2 * V)
    Y.list[[i]] <- X %*% theta_true + epsilon
  }

  # ---------------------------------------------------------------------------
  # Return
  # ---------------------------------------------------------------------------
  return(list(
    X               = X,
    Y.list          = Y.list,
    optimal_x       = optimal_x,
    optimal_weights = optimal_weights
  ))
}

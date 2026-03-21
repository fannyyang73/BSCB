compute_NG_param <- function(X, Y, V, hyperparameter = c("empirical", "unit_info", "g_prior")) {

  hyperparameter <- match.arg(hyperparameter)

  # ── 1. dataset ──────────────────────────────────────────

  n <- nrow(X)
  q <- ncol(X)          # q = p + 1 (intercept included)

  # ── 2. least squares estimator ──────────────────────────────
  V_inv    <- solve(V)
  XtVX     <- t(X) %*% V_inv %*% X        # (q x q)
  XtVY     <- t(X) %*% V_inv %*% Y
  XtVX_inv <- solve(XtVX)

  beta_hat <- XtVX_inv %*% XtVY  # OLS

  residual  <- Y - X %*% beta_hat
  RSS       <- as.numeric(t(residual) %*% V_inv %*% residual)
  sigma2_hat <- RSS / (n - q)               # OLS of variance

  # ── 3. determint the hyperparameters ───────────────────────
  if (hyperparameter == "empirical") {
    # Normal-Gamma with empirical hyperparameters
    # P = c * X'V^{-1}X,  mu_0 = beta_hat
    c_value <- 0.001
    mu_0    <- beta_hat
    cal_P   <- c_value * diag(q)          # prior precision matrix
    alpha_0 <- 1
    beta_0  <- sigma2_hat

  } else if (hyperparameter == "unit_info") {
    # Unit Information Prior
    # P = (1/n) * X'V^{-1}X,
    mu_0    <- beta_hat
    cal_P   <- (1 / n) * XtVX
    alpha_0 <- 0.5
    beta_0  <-  sigma2_hat / 2

  } else if (hyperparameter == "g_prior") {
    # Zellner's g-Prior, g = n
    g       <- n
    mu_0    <-  beta_hat #matrix(0, nrow = q, ncol = 1) #beta_hat
    cal_P   <- (1 / g) * XtVX            # P = (1/g) * X'V^{-1}X
    alpha_0 <- 0.5
    beta_0  <- sigma2_hat / 2

  }
  # ── 4. Compute the updated posterior parameters ───────────

  inv_cal_P      <- solve(cal_P)

  # (1) conditional posterior of theta | tau, Y
  # theta | tau, Y ~ Normal(mu_n, tau^{-1}*inv_cal_Pn)
  cal_Pn <- XtVX + cal_P
  inv_cal_Pn <- solve(cal_Pn)
  mu_n <- inv_cal_Pn %*% ( XtVY + cal_P %*% mu_0)

  # (2) marginal posterior of tau | Y
  # tau | Y ~ Gamma(alpha_n, beta_n)
  alpha_n <- (n+2*alpha_0)/2
  c_scalar <- t(Y) %*% V_inv %*% Y + t(mu_0) %*% cal_P %*%  mu_0 - t(mu_n) %*% cal_Pn %*% mu_n
  beta_n <- beta_0 + (c_scalar/2)

  # (3) marginal posterior of theta | Y
  # theta | Y ~ multivariate-t(dof_NG, mu_NG, scale_NG )
  dof_NG <- 2*alpha_n
  mu_NG  <- mu_n
  scale_NG <- inv_cal_Pn * as.numeric(beta_n/alpha_n)

  # (4) prior predictive of Y
  # Y ~ multivariate-t(dof_m,  mu_m,  scale_m)
  dof_m   <- 2 * alpha_0
  mu_m  <- X %*% mu_0                                  # (n x 1)
  scale_m <- (beta_0 / alpha_0) * (V + X %*% inv_cal_P %*% t(X))  # (n x n)

  # (5)  posterior predictive of Y
  # Y_rep ~ multivariate-t(dof_p, mu_p, scale_p)
  dof_p <- 2*alpha_n
  mu_p <- X %*% mu_n
  H <- X %*% inv_cal_Pn %*% t(X) + V
  scale_p <- H * as.numeric(beta_n/alpha_n)




  # ── 5. Return ────────────────────────────────────────────
  return(list(
    hyperparameter = hyperparameter,
    # initial hyper parameter
    mu_0           = mu_0,
    cal_P          = cal_P,
    alpha_0        = alpha_0,
    beta_0         = beta_0,
    # updated hyper parameter
    mu_n           = mu_n,
    cal_Pn          = cal_Pn,
    alpha_n        = alpha_n,
    beta_n         = beta_n,

    # marginal posterior of theta
    marginal_pos_theta_dof = dof_NG,
    marginal_pos_theta_mean = mu_NG,
    marginal_pos_theta_scale = scale_NG,

    # posterior predictive of Y
    pos_pred_Y_dof  = dof_p,
    pos_pred_Y_mean  = mu_p,
    pos_pred_Y_scale  = scale_p,


    # prior predictive/ marginal predictive of Y
    prior_pred_Y_dof   = dof_m,
    prior_pred_Y_mean  = mu_m,
    prior_pred_Y_scale = scale_m,

    # OLS
    beta_hat       = beta_hat,
    sigma2_hat     = sigma2_hat,
    RSS            = RSS,
    n              = n,
    q              = q
  ))
}

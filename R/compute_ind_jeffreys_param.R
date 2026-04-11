compute_IJ_param <- function(X, Y, V) {


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
  RSS       <- as.numeric(t(residual) %*% V_inv %*% residual) # SSE
  sigma2_hat <- RSS / (n - q)               # OLS of variance

  # ── 3. Compute the updated posterior parameters ───────────

  # (3) marginal posterior of theta | Y
  # theta | Y ~ multivariate-t(dof_NG, mu_NG, scale_NG )
  dof_J <- n - q
  mu_J  <- beta_hat
  scale_J <- sigma2_hat * XtVX_inv



  # ── 4. Return ────────────────────────────────────────────
  return(list(

    # marginal posterior of theta
    marginal_pos_theta_dof = dof_J,
    marginal_pos_theta_mean = mu_J,
    marginal_pos_theta_scale = scale_J,


    # OLS
    beta_hat       = beta_hat,
    sigma2_hat     = sigma2_hat,
    RSS            = RSS,
    n              = n,
    q              = q
  ))
}

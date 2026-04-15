data {
  int<lower=1> N;
  int<lower=2> p;          // polynomial degree (2 or 3)
  vector[N] x;
  vector[N] y;
  matrix[N, N] V;
  real<lower=0> theta_prior_sd;
  vector[p+1] theta_prior_mean;
  int<lower=0, upper=1> use_likelihood;
  int<lower=0, upper=1> HMC_prior_type;
  real<lower=0> normal_theta_sd;
  real<lower=0> normal_sigma_sd;
  real<lower=0> cauchy_scale;
}
transformed data {
  cholesky_factor_cov[N] L_V;
  matrix[N, p+1] X_mat;
  L_V = cholesky_decompose(V);
  for (i in 1:N) {
    for (k in 0:p) {
      X_mat[i, k+1] = x[i]^k;
    }
  }
}
parameters {
  vector[p+1] theta;
  real<lower=0> sigma;
}
transformed parameters {
  vector[N] mu;
  mu = X_mat * theta;
}
model {
  if (HMC_prior_type == 0) {
    theta ~ normal(0, normal_theta_sd);
    sigma ~ normal(0, normal_sigma_sd);
  } else {
    theta ~ normal(theta_prior_mean, theta_prior_sd);
    sigma ~ cauchy(0, cauchy_scale);
  }
  if (use_likelihood == 1)
    y ~ multi_normal_cholesky(mu, sigma * L_V);
}
generated quantities {
  vector[N] Y_pred;
  Y_pred = multi_normal_cholesky_rng(mu, sigma * L_V);
}

data {
  int<lower=0> N;
  vector[N] x;
  vector[N] y;
  matrix[N, N] V;
  real<lower=0> theta_prior_sd; 
  real theta_prior_mean0;
  real theta_prior_mean1;
  real theta_prior_mean2;
  int<lower=0, upper=1> use_likelihood;  // 0 = prior predictive, 1 = posterior
  int<lower=0, upper=1> HMC_prior_type;
  real<lower=0> normal_theta_sd; 
  real<lower=0> normal_sigma_sd; 
  real<lower=0> cauchy_scale;
}
transformed data {
  cholesky_factor_cov[N] L_V;
  matrix[N, 3] X_mat;
  L_V = cholesky_decompose(V);
  for (i in 1:N) {
    X_mat[i, 1] = 1.0;
    X_mat[i, 2] = x[i];
    X_mat[i, 3] = square(x[i]);
  }
}
parameters {
  vector[3] theta;       
  real<lower=0> sigma;
}
transformed parameters {
  vector[N] mu;
  mu = X_mat * theta;    
}
model {
  if (HMC_prior_type == 0){ // 0: normal-normal
    theta[1] ~ normal(0, normal_theta_sd);
    theta[2] ~ normal(0, normal_theta_sd);
    theta[3] ~ normal(0, normal_theta_sd);
    sigma  ~ normal(0, normal_sigma_sd);
  }
  else if (HMC_prior_type == 1){ // 1: normal-halfCauchy
    theta[1] ~ normal(theta_prior_mean0, theta_prior_sd);
    theta[2] ~ normal(theta_prior_mean1, theta_prior_sd);
    theta[3] ~ normal(theta_prior_mean2, theta_prior_sd);
    sigma  ~ cauchy(0, cauchy_scale);
  }
  
  
  if (use_likelihood == 1)
    y ~ multi_normal_cholesky(mu, sigma * L_V);
}
generated quantities {
  vector[N] Y_pred;
  Y_pred = multi_normal_cholesky_rng(mu, sigma * L_V);
}


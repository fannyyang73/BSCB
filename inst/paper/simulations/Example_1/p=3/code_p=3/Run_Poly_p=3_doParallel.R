######################################################################
### Bayesian Simultaneous Credible Bands for Polynomial Regression ###
######################################################################
## Author: Fei Yang ##
## 2026-03-24
## lambda by find_global_maximum()
## Hyperparameters have been adjusted
## Re-proofed version
## Adjust x_scale to only mean-centered
## Test empirical coverage rate for BSCB
## Fix X across all batches
## Add D_optim
## Add 2 ways: coverage probability & empirical coverage rate
## Foreach, inside the loop use HMC
## Add the Jeffery's non-informative prior
## dof <- n
## read in the datasets of the HMC
## Order is 2
## Errors could be iid or with auto-regressive structure
## Add set.seed
## Add independent Jefferys prior
## Add the unit-information prior and the g-prior for the normal-gamma prior
## change HMC from normal-normal to normal-half Cauchy(0,2)
## HPC pool
## h_optim
## doParallel
## add BPCB-NC


####################################.     Read in the datasets   ####################################
load("Data_CenterStatus1_DesignMat2_Batch1.RData")
source("Functions_Poly_p=3.R")

library(doParallel)
library(foreach)
# library(doMPI)
# library(Rmpi)
library(psych)
library(mvtnorm) # For multivariate normal/t distribution
library(EnvStats)
#library(caret)
library(MASS)
library(glmnet)
library(OptimalDesign)
library(Metrics)
library(combinat)
library(AlgDesign)
library(cmdstanr)
library(posterior)


theta <- c(1, 2, -1, 0.5)
theta_true <- c(1, 2, -1, 0.5)
rho <- 0.1
AR_index <- 0 # 0: iid error; 1: auto-regressive error

p <- 3
e_sd <- 0.2
n <- 20

replication <- 1000
run_once <- 0 # 1: only run once for the picked id; 0: run for all replications
pick_id <- 37

L <- 500000
alpha <- 1 #values due to non-informative prior
n_new <- 1000
a <- -5
b <- 5
delta <- 0.05
c_value <- 0.001
draw_num <- 10000

if(AR_index == "0"){ # 0: iid error; 1: auto-regressive error
  V <- diag(n)
}else if(AR_index == "1"){V <- rho^abs(outer(1:n, 1:n, "-")) / (1 - rho^2)}
V_inv <- solve(V)

HMC_setting <- 0 # 0: Use the normal-gamma prior; 1:Use HMC; 2: Use the non-informative prior; 3: Use the independent Jefferys prior
hyper_prior_setting <- 2 # 0: empirical ; 1: unit-information; 2: g-prior
HMC_prior_type <- 1 # 0: normal-normal; 1: normal-halfCauchy
optim_type <- 1 # 0: gloabl max; 1: Liu's method
normal_theta_sd <- 10
normal_sigma_sd  <- 5
cauchy_scale <-  2
mu_setting <- 1 # 1:OLS; 2:Huber; 3:Ridge
P_setting <- 1 #1:I_p; 2: eta
center_index <- 1 # 1:centered; 0:not centered; 2: standardization
design_index <- 2 # 1:Equal space; 2:D-optimal; 3: D-optimal with equal weights and pre-specified levels
x_setting <- 1 # 1:uniform distribution; 2:standard normal distribution


##### design_index = 1 (ES) #####
n_ES_x <- 20
##### design_index = 2 (DO) #####
n_DO_init_x <- 300000 # large values guarantees 6 levels in produced x
##### design_index = 3 (DOE) #####
n_DOE_init_x <- 20 # initial levels in x; n_DOE_init_x should be larger than x_level
x_level <- 10 # final levels in x; Must make sure that n/x_level should be an integer
##### design_index = 4 (DOS) #####
n_k <- 13 # Must be an odd number
tail_weight <- 1.8

N_index <- 1 # Batch_index
num_cores <- 6

########################       Parallel Computing       ########################
multiResultClass <- function(i=NULL, lambda_best_optim=NULL, coverage_flag=NULL, mu_star=NULL, cov_theta=NULL,  coverage_draw=NULL, seed_used=NULL, elapsed_time=NULL, theta_quantiles=NULL, theta_mat=NULL){
  me <- list(i = i, lambda_best_optim = lambda_best_optim, coverage_flag = coverage_flag, mu_star = mu_star,
             cov_theta=cov_theta,  coverage_draw = coverage_draw, seed_used = seed_used, elapsed_time = elapsed_time, theta_quantiles = theta_quantiles, theta_mat=theta_mat )
  class(me) <- append(class(me),"multiResultClass") ## Set the name for the class
  return(me)
}


set.seed(2025)
replication_seeds <- sample.int(.Machine$integer.max, replication)
cl <- makeCluster(num_cores)
registerDoParallel(cl)
# Export required functions & objects to worker nodes
clusterExport(cl, c("X", "Y.list", "n", "p", "L", "a", "b", "delta", "alpha",
                    "tr", "rmvt", "f_L_SCB", "f_U_SCB", "theta_true", "L_SCB", "U_SCB", "prob_x_SCB", "truncated_normal_density",
                    "normal_density", "uniform_density", "integrand_SCB",
                    "prob_x_SCB_HMC", "integrand_SCB_HMC","find_global_maximum","fn_Bayes_PCP", "fn_Bayes_ECR", "fn_Freq_ECR",
                    "find_global_maximum_h_all","sup_T_Bayes_PSCP","sup_T_Bayes_ESCR","sup_T_simFSCB","sup_T_Freq_ESCR",
                    "replication_seeds","compute_NG_param","cauchy_scale","HMC_prior_type","normal_theta_sd","normal_sigma_sd"))


start_time <- Sys.time()
result_mat <- matrix(NA, nrow = replication, ncol = 6)
mu_star_mat <- matrix(nrow = replication, ncol = p+1)
cov_theta_list <- list()
coverage_prob_list <- list()
lambda_best_optim_list <- list()
coverage_draw_list <- list()
seed_used_list <- list()
elapsed_time_list <- list()
theta_quantiles_list <- list()
theta_mat_list <- list()



if(run_once == 1){
  iter_range <- pick_id   # 只跑1次
} else {
  iter_range <- 1:replication
}

results <- foreach(i = iter_range,  .packages = c("mvtnorm","MASS","glmnet","cmdstanr", "posterior")) %dopar% { #, "caret"
  cmdstanr::set_cmdstan_path()
  rep_start_time <- Sys.time()

  set.seed(replication_seeds[i]) # set.seed to ensure reproducibility
  Y <- Y.list[[i]]
  lambda_sup <- numeric()
  lambda_draw <- numeric()

if(HMC_setting == "1"){
  ##########################################################################################
  ############.                             (1) Use HMC.                        ############
  ##########################################################################################

  worker_mod <- cmdstan_model("HMC_p=3.stan", compile = TRUE, force_recompile = TRUE)

  dat_for_use <- list(X = X, V = V, Y = Y.list)
  E_NG <- compute_NG_param(data = dat_for_use, pick_id = i, prior = "empirical")

  stan_data <- list(N = n, x = X[,2], y = as.numeric(Y.list[[i]]),
                    V = V,
                    theta_prior_sd    = sqrt(E_NG$beta_0),
                    theta_prior_mean0 = E_NG$beta_hat[1,],
                    theta_prior_mean1 = E_NG$beta_hat[2,],
                    theta_prior_mean2 = E_NG$beta_hat[3,],
                    theta_prior_mean3 = E_NG$beta_hat[4,],
                    use_likelihood    = 1,
                    HMC_prior_type    = HMC_prior_type,
                    normal_theta_sd   = normal_theta_sd,
                    normal_sigma_sd   = normal_sigma_sd,
                    cauchy_scale      = cauchy_scale)
  fit <- worker_mod$sample(
    data = stan_data,
    iter_sampling = 4000,
    iter_warmup = 4000,
    chains = 4,
    seed = 123,
    parallel_chains = 1,
    refresh = 0,
    adapt_delta = 0.95,
    max_treedepth = 15,
    save_cmdstan_config = TRUE
  )

  # Extract draws
  posterior_samples <- fit$draws(format = "df")
  # Sample size
  n_samples <- nrow(posterior_samples)
  # Extract and organize theta samples
  theta_mat <- cbind(
    posterior_samples$`theta[1]`,
    posterior_samples$`theta[2]`,
    posterior_samples$`theta[3]`,
    posterior_samples$`theta[4]`
  )

  # Store results
  mu_star <- colMeans(theta_mat)
  mu_star_mat[i, ] <- mu_star
  cov_theta <- cov(theta_mat)
  cov_theta_list[[i]] <- cov_theta

  theta_mat_list[[i]] <- theta_mat
  theta_quantiles <- apply(theta_mat, 2, quantile, probs = c(delta/2, 1 - delta/2))
  rownames(theta_quantiles) <- c(paste0("Q", delta/2*100, "%"), paste0("Q", (1-delta/2)*100, "%"))
  colnames(theta_quantiles) <- c("theta1", "theta2", "theta3","theta4")


  #####################.  Compute the critical constant  #####################
  for (j in 1:L) {
    iii <- sample(1:n_samples, 1)
    theta_hat <- theta_mat[iii,]

    if(optim_type == "0"){
      result <- find_global_maximum(fn = fn_Bayes_PCP, a, b, order_form, theta = theta_hat, mu_star = mu_star, cov_mat = cov_theta)
    }else if (optim_type == "1"){
      result <- sup_T_Bayes_PSCP(a, b, theta_hat, mu_star, cov_mat = cov_theta)
    }
    lambda_sup[j] <- result$maximum
  }

  lambda_best_optim <- quantile(lambda_sup, probs = 1 - delta)
  cat("The best lambda for the", i, "round is", lambda_best_optim, "\n")

  #####################.  Empirical Coverage Rate.  #####################


  if(optim_type == "0"){
    result_ECR <- find_global_maximum(fn = fn_Bayes_ECR, a, b, order_form, theta = theta_true, mu_star = mu_star, cov_mat = cov_theta)
  }else if (optim_type == "1"){
    result_ECR <- sup_T_Bayes_ESCR(a, b, theta_true, mu_star, cov_mat = cov_theta)
  }
  coverage_flag <- ifelse( lambda_best_optim >= result_ECR$maximum, 1, 0)
  if(coverage_flag == 0){cat("Round",i,"is rejected by the empirical coverage rate\n")}


  #####################.  Bayesian Coverage Probability  #####################

  cover_num <- 0
  for (j in 1:draw_num) {
    S <- nrow(theta_mat)
    idx <- sample(1:S, 1)
    theta_hat <- theta_mat[idx, ]

    if(optim_type == "0"){
      result_draw <- find_global_maximum(fn = fn_Bayes_PCP, a, b, order_form, theta = theta_hat, mu_star = mu_star, cov_mat = cov_theta)
    }else if (optim_type == "1"){
      result_draw <- sup_T_Bayes_PSCP(a, b, theta_hat, mu_star, cov_mat = cov_theta)
    }

    lambda_draw[j] <- result_draw$maximum

    if(lambda_draw[j] <= lambda_best_optim){
      cover_num <- cover_num + 1
    }
  }
  coverage_draw <- cover_num/draw_num
  cat("The Bayesian coverage rate of Bayesian simultaneous credible band for Round", i , "is:", coverage_draw, "\n")



}else if(HMC_setting == "0"){
  ##########################################################################################
  ############.         (0) Use the normal- gamma conjugate prior               ############
  ##########################################################################################

  data <- list(X = X, V = V, Y = Y.list)

  if(hyper_prior_setting == "0"){ # 0: "empirical", 1: "unit_info", 2:"g_prior"
    E_NG <- compute_NG_param(data = data, pick_id = i, prior="empirical")
    mu_star <-   E_NG$marginal_pos_theta_mean
    dof <- E_NG$marginal_pos_theta_dof
    scale_mat <- E_NG$marginal_pos_theta_scale
  }else if (hyper_prior_setting == "1"){
    U_NG <- compute_NG_param(data = data, pick_id = i, prior="unit_info")
    mu_star <-   U_NG$marginal_pos_theta_mean
    dof <- U_NG$marginal_pos_theta_dof
    scale_mat <- U_NG$marginal_pos_theta_scale
  }else if (hyper_prior_setting == "2"){
    G_NG <- compute_NG_param(data = data, pick_id = i, prior="g_prior")
    mu_star <-   G_NG$marginal_pos_theta_mean
    dof <- G_NG$marginal_pos_theta_dof
    scale_mat <- G_NG$marginal_pos_theta_scale
  }
  the_value <- (dof) / (dof - 2)
  cov_theta <- drop(the_value) * scale_mat

  #####################.  Compute the critical constant  #####################
  for (j in 1:L) {
    theta_hat <- rmvt(n = 1, delta = mu_star, sigma = scale_mat, df = dof, type = "shifted")

    if(optim_type == "0"){
      result <- find_global_maximum(fn = fn_Bayes_PCP, a, b, order_form, theta = theta_hat, mu_star = mu_star, cov_mat = cov_theta)
    }else if (optim_type == "1"){
      result <- sup_T_Bayes_PSCP(a, b, theta_hat, mu_star, cov_mat = cov_theta)
    }
    lambda_sup[j] <- result$maximum
  }


  lambda_best_optim <- quantile(lambda_sup, probs = 1 - delta)
  cat("The best lambda for the", i, "round is",lambda_best_optim, "\n")

  #####################.  Empirical Coverage Rate.  #####################

  if(optim_type == "0"){
    result_ECR <- find_global_maximum(fn = fn_Bayes_ECR, a, b, order_form, theta = theta_true, mu_star = mu_star, cov_mat = cov_theta)
  }else if (optim_type == "1"){
    result_ECR <- sup_T_Bayes_ESCR(a, b, theta_true, mu_star, cov_mat = cov_theta)
  }
  coverage_flag <- ifelse( lambda_best_optim >= result_ECR$maximum, 1, 0)
  if(coverage_flag == 0){cat("Round",i,"is rejected by the empirical coverage rate\n")}



  #####################.  Bayesian Coverage Probability  #####################
  cover_num <- 0
  for (j in 1:draw_num) {
    theta_hat <- rmvt(n = 1, delta = mu_star, sigma = scale_mat, df = dof, type = "shifted")

    if(optim_type == "0"){
      result_draw <- find_global_maximum(fn = fn_Bayes_PCP, a, b, order_form, theta = theta_hat, mu_star = mu_star, cov_mat = cov_theta)
    }else if (optim_type == "1"){
      result_draw <- sup_T_Bayes_PSCP(a, b, theta_hat, mu_star, cov_mat = cov_theta)
    }
    lambda_draw[j] <- result_draw$maximum

    if(lambda_draw[j] <= lambda_best_optim){
      cover_num <- cover_num + 1
    }
  }
  coverage_draw <- cover_num/draw_num
  cat("The Bayesian coverage rate of Bayesian simultaneous credible band for Round", i , "is:", coverage_draw, "\n")


}else if(HMC_setting == "2"){
  ##########################################################################################
  ############.       (2) Use the Jeffery's non-informative prior.              ############
  ##########################################################################################
  XtX <- t(X) %*% X
  XtX_inv <- solve(XtX)
  XtVX <- t(X) %*% V_inv %*% X

  mu_star <- solve(t(X)%*%V_inv%*%X)%*%t(X)%*%V_inv%*%Y
  error <- Y - X %*% mu_star

  variance <- t(error) %*% V_inv %*% error / n

  dof <- n

  scale_mat <- drop(variance) * solve(XtVX)

  the_value <- (dof) / (dof - 2)
  cov_theta <- drop(the_value) * scale_mat

  #####################.  Compute the critical constant  #####################
  for (j in 1:L) {
    theta_hat <- rmvt(n = 1, delta = mu_star, sigma = scale_mat, df = dof, type = "shifted")

    if(optim_type == "0"){
      result <- find_global_maximum(fn = fn_Bayes_PCP, a, b, order_form, theta = theta_hat, mu_star = mu_star, cov_mat = cov_theta)
    }else if (optim_type == "1"){
      result <- sup_T_Bayes_PSCP(a, b, theta_hat, mu_star, cov_mat = cov_theta)
    }
    lambda_sup[j] <- result$maximum

  }

  lambda_best_optim <- quantile(lambda_sup, probs = 1 - delta)
  cat("The best lambda for the", i, "round is",lambda_best_optim, "\n")

  #####################.  Empirical Coverage Rate.  #####################

  if(optim_type == "0"){
    result_ECR <- find_global_maximum(fn = fn_Bayes_ECR, a, b, order_form, theta = theta_true, mu_star = mu_star, cov_mat = cov_theta)
  }else if (optim_type == "1"){
    result_ECR <- sup_T_Bayes_ESCR(a, b, theta_true, mu_star, cov_mat = cov_theta)
  }
  coverage_flag <- ifelse( lambda_best_optim >= result_ECR$maximum, 1, 0)
  if(coverage_flag == 0){cat("Round",i,"is rejected by the empirical coverage rate\n")}




  #####################.  Bayesian Coverage Probability  #####################
  cover_num <- 0
  for (j in 1:draw_num) {
    theta_hat <- rmvt(n = 1, delta = mu_star, sigma = scale_mat, df = dof, type = "shifted")

    if(optim_type == "0"){
      result_draw <- find_global_maximum(fn = fn_Bayes_PCP, a, b, order_form, theta = theta_hat, mu_star = mu_star, cov_mat = cov_theta)
    }else if (optim_type == "1"){
      result_draw <- sup_T_Bayes_PSCP(a, b, theta_hat, mu_star, cov_mat = cov_theta)
    }
    lambda_draw[j] <- result_draw$maximum

    if(lambda_draw[j] <= lambda_best_optim){
      cover_num <- cover_num + 1
    }
  }
  coverage_draw <- cover_num/draw_num
  cat("The Bayesian coverage rate of Bayesian simultaneous credible band for Round", i , "is:", coverage_draw, "\n")



}else if(HMC_setting == "3"){
  ##########################################################################################
  #########.       (3) Use the independent Jeffery's non-informative prior.        #########
  ##########################################################################################
  XtX <- t(X) %*% X
  XtX_inv <- solve(XtX)
  XtVX <- t(X) %*% V_inv %*% X

  mu_star <- solve(t(X)%*%V_inv%*%X)%*%t(X)%*%V_inv%*%Y
  error <- Y - X %*% mu_star
  dof <- n - p - 1
  variance <- t(error) %*% V_inv %*% error / dof

  scale_mat <- drop(variance) * solve(XtVX)

  the_value <- (dof) / (dof - 2)
  cov_theta <- drop(the_value) * scale_mat

  #####################.  Compute the critical constant  #####################
  for (j in 1:L) {
    theta_hat <- rmvt(n = 1, delta = mu_star, sigma = scale_mat, df = dof, type = "shifted")

    if(optim_type == "0"){
      result <- find_global_maximum(fn = fn_Bayes_PCP, a, b, order_form, theta = theta_hat, mu_star = mu_star, cov_mat = cov_theta)
    }else if (optim_type == "1"){
      result <- sup_T_Bayes_PSCP(a, b, theta_hat, mu_star, cov_mat = cov_theta)
    }
    lambda_sup[j] <- result$maximum

  }

  lambda_best_optim <- quantile(lambda_sup, probs = 1 - delta)
  cat("The best lambda for the", i, "round is",lambda_best_optim, "\n")

  #####################.  Empirical Coverage Rate.  #####################

  if(optim_type == "0"){
    result_ECR <- find_global_maximum(fn = fn_Bayes_ECR, a, b, order_form, theta = theta_true, mu_star = mu_star, cov_mat = cov_theta)
  }else if (optim_type == "1"){
    result_ECR <- sup_T_Bayes_ESCR(a, b, theta_true, mu_star, cov_mat = cov_theta)
  }
  coverage_flag <- ifelse( lambda_best_optim >= result_ECR$maximum, 1, 0)
  if(coverage_flag == 0){cat("Round",i,"is rejected by the empirical coverage rate\n")}




  #####################.  Bayesian Coverage Probability  #####################
  cover_num <- 0
  for (j in 1:draw_num) {
    theta_hat <- rmvt(n = 1, delta = mu_star, sigma = scale_mat, df = dof, type = "shifted")

    if(optim_type == "0"){
      result_draw <- find_global_maximum(fn = fn_Bayes_PCP, a, b, order_form, theta = theta_hat, mu_star = mu_star, cov_mat = cov_theta)
    }else if (optim_type == "1"){
      result_draw <- sup_T_Bayes_PSCP(a, b, theta_hat, mu_star, cov_mat = cov_theta)
    }
    lambda_draw[j] <- result_draw$maximum

    if(lambda_draw[j] <= lambda_best_optim){
      cover_num <- cover_num + 1
    }
  }
  coverage_draw <- cover_num/draw_num
  cat("The Bayesian coverage rate of Bayesian simultaneous credible band for Round", i , "is:", coverage_draw, "\n")

}

  rep_end_time <- Sys.time()
  elapsed_time <- as.numeric(difftime(rep_end_time, rep_start_time, units = "secs"))


result <- multiResultClass()
result$i <- i
result$lambda_best_optim <- lambda_best_optim
result$coverage_flag <- coverage_flag
result$mu_star <- mu_star
result$cov_theta <- cov_theta
result$coverage_draw <- coverage_draw
result$seed_used <- replication_seeds[i]
result$elapsed_time <- elapsed_time
if(HMC_setting == "1"){
  result$theta_quantiles <- theta_quantiles
  result$theta_mat <- theta_mat
}
return(result)
}

stopCluster(cl)
end_time <- Sys.time()
duration_cov <- end_time - start_time


# if(HMC_setting == 1){
#   job_outdir <- file.path(
#     "/mnt/iusers01/maths01/w85950fy/scratch",
#     Sys.getenv("SLURM_JOB_ID", "local")
#   )
#   if(dir.exists(job_outdir)){
#     unlink(job_outdir, recursive = TRUE)
#     cat("Cleaned up Stan output directory:", job_outdir, "\n")
#   }
# }

if(run_once == 1){
  collect_range <- 1  # results 只有1个元素
} else {
  collect_range <- 1:replication
}

# for(i in 1:replication){
#   result_mat[i,1] <- results[[i]]$i
#   result_mat[i,2] <- results[[i]]$lambda_best_optim
#   result_mat[i,3] <- results[[i]]$coverage_flag
#   result_mat[i,4] <- results[[i]]$coverage_draw
#   result_mat[i,5] <- results[[i]]$seed_used
#   result_mat[i,6] <- results[[i]]$elapsed_time
#
#   mu_star_mat[i,] <- results[[i]]$mu_star
#   cov_theta_list[[i]] <- results[[i]]$cov_theta
#   lambda_best_optim_list[[i]] <- results[[i]]$lambda_best_optim
#   coverage_draw_list[[i]] <- results[[i]]$coverage_draw
#   seed_used_list[[i]] <- results[[i]]$seed_used
#   elapsed_time_list[[i]] <- results[[i]]$elapsed_time
#   if(HMC_setting == "1"){
#     theta_quantiles_list[[i]] <- results[[i]]$theta_quantiles
#     theta_mat_list[[i]] <- results[[i]]$theta_mat}
# }
for(k in collect_range){
  store_id <- if(run_once == 1) pick_id else k

  result_mat[store_id, 1] <- results[[k]]$i
  result_mat[store_id, 2] <- results[[k]]$lambda_best_optim
  result_mat[store_id, 3] <- results[[k]]$coverage_flag
  result_mat[store_id, 4] <- results[[k]]$coverage_draw
  result_mat[store_id, 5] <- results[[k]]$seed_used
  result_mat[store_id, 6] <- results[[k]]$elapsed_time
  mu_star_mat[store_id, ] <- results[[k]]$mu_star
  cov_theta_list[[store_id]] <- results[[k]]$cov_theta
  lambda_best_optim_list[[store_id]] <- results[[k]]$lambda_best_optim
  if(HMC_setting == 1){
    theta_quantiles_list[[store_id]] <- results[[k]]$theta_quantiles
    theta_mat_list[[store_id]] <- results[[k]]$theta_mat
  }
}

colnames(result_mat) <- c("Index", "Lambda_best_optim", "Empirical_Coverage_Rate", "Coverage_Rate_Draw","Seed_used","Elapsed_time (seconds)")



# Save results
write.csv(result_mat, "current_result.csv", row.names = FALSE)

if(run_once == 1){
  sucess_count_SCB <- result_mat[pick_id, 3]
  SCB_emp_cov_rate <- sucess_count_SCB
  SCB_avg_cov_draw <- result_mat[pick_id, 4]
  MSE <- sum((theta_true - mu_star_mat[pick_id, ])^2)
} else {
  sucess_count_SCB <- sum(result_mat[, 3], na.rm = TRUE)
  SCB_emp_cov_rate <- sucess_count_SCB / replication
  SCB_avg_cov_draw <- mean(result_mat[, 4], na.rm = TRUE)
  ac_mat <- t(replicate(replication, theta_true))
  MSE <- mean(rowSums((ac_mat - mu_star_mat)^2, na.rm = TRUE), na.rm = TRUE)
}

# sucess_count_SCB <- sum(result_mat[, 3]) # Count success cases
# SCB_emp_cov_rate <- sucess_count_SCB / replication  # Empirical Coverage Rate
# SCB_avg_cov_draw <- mean(result_mat[,4])  # Posterior Coverage Rate
# ac_mat <- t(replicate(replication, theta_true))
# MSE <- mean(rowSums((ac_mat-mu_star_mat)^2))

final_result <- matrix(nrow=1, ncol = 5)
final_result[1,1] <- N_index
final_result[1,2] <- SCB_emp_cov_rate
final_result[1,3] <- SCB_avg_cov_draw
final_result[1,4] <- MSE
final_result[1,5] <- duration_cov
colnames(final_result) <- c("Batch_index", "SCB_emp_cov_rate", "SCB_avg_cov_draw","MSE", "Total Time")
write.csv(final_result, paste0("final_result_batch",N_index,".csv"), row.names = FALSE )

cat("The Empirical Coverage Rate for the Bayesian simultaneous credible band is:", SCB_emp_cov_rate, "\n")
cat("The Average Coverage Rate (by drawing) for the Bayesian simultaneous credible band is:", SCB_avg_cov_draw, "\n")
cat("The MSE for theta in total of",replication, "rounds is:", MSE, "\n")
print(duration_cov)


save.image(file = "All_result.RData")
# mpi.quit()


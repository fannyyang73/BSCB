##########################################################################################
####################################.     Functions   ####################################
##########################################################################################

#' Create polynomial basis function
#' @keywords internal
#' @export
create_order_form <- function(p) {
  if (!is.numeric(p) || p < 0 || p != floor(p)) {
    stop("p must be a non-negative integer")
  }

  if (p == 1) {
    function(x) {
      if (is.matrix(x) || length(x) > 1) {
        return(cbind(1, x))
      } else {
        return(c(1, x))
      }
    }
  } else if (p == 2) {
    function(x) {
      if (is.matrix(x) || length(x) > 1) {
        return(cbind(1, x, x^2))
      } else {
        return(c(1, x, x^2))
      }
    }
  } else if (p == 3) {
    function(x) {
      if (is.matrix(x) || length(x) > 1) {
        return(cbind(1, x, x^2, x^3))
      } else {
        return(c(1, x, x^2, x^3))
      }
    }
  } else {
    # 通用版本
    powers <- 0:p
    function(x) {
      result <- outer(x, powers, `^`)
      if (length(x) == 1) return(as.vector(result))
      return(result)
    }
  }
}


#########################################################
#########.             Optimization           ###########
#########################################################

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

fn_Bayes_ECR <- function(x, theta_true, mu_star, cov_theta) {
  x_i <- order_form(x)
  numerator <- abs((x_i)%*%t(theta_true-t(mu_star)))
  denominator <- sqrt(x_i%*%cov_theta%*%t(t(x_i)))
  lambda <- numerator%*%solve(denominator)
  return(lambda)
}


fn_Freq_ECR <- function(x, theta_true, lm_theta_hat, S, inv) {
  x_i <- order_form(x)
  numerator <- abs((x_i) %*% t(theta_true - t(lm_theta_hat)))
  denominator <- sqrt(x_i%*%(as.numeric(S^2) * inv)%*%t(t(x_i))) # S, inv
  lambda <- numerator %*% solve(denominator)
  return(lambda)
}

#### find the global maximum of T(x) form function
find_global_maximum <- function(fn, a, b, order_form, theta, mu_star, cov_mat,
                                tol = 1e-6, n_grid = 100) {

  # 定义分子函数（不取绝对值）
  numerator_fn <- function(x) {
    x_i <- order_form(x)
    as.numeric((x_i) %*% t(theta - t(mu_star)))
  }

  # 定义分母函数
  denominator_fn <- function(x) {
    x_i <- order_form(x)
    as.numeric(sqrt(x_i %*% cov_mat %*% t(t(x_i))))
  }

  # 存储候选点
  candidate_x <- c()
  candidate_values <- c()

  # (1) 两个端点
  candidate_x <- c(candidate_x, a, b)
  candidate_values <- c(candidate_values, fn(a), fn(b))

  # (2) 寻找分子为零的点（尖点）
  # 使用网格搜索 + uniroot
  grid <- seq(a, b, length.out = n_grid)
  num_values <- sapply(grid, numerator_fn)

  # 检测符号变化
  sign_changes <- which(diff(sign(num_values)) != 0)

  for (idx in sign_changes) {
    tryCatch({
      zero_point <- uniroot(numerator_fn,
                            interval = c(grid[idx], grid[idx + 1]),
                            tol = tol)$root
      candidate_x <- c(candidate_x, zero_point)
      candidate_values <- c(candidate_values, fn(zero_point))
    }, error = function(e) {})
  }

  # (3) 寻找导数为零的点
  # 使用数值导数
  fn_deriv <- function(x, h = 1e-6) {
    (fn(x + h) - fn(x - h)) / (2 * h)
  }

  # 在每个子区间寻找导数为零的点
  # 排除已找到的尖点附近
  search_intervals <- list()

  # 将[a,b]分割，避开尖点
  zero_points <- candidate_x[candidate_x > a & candidate_x < b]
  zero_points <- sort(zero_points)

  if (length(zero_points) == 0) {
    search_intervals <- list(c(a, b))
  } else {
    # 创建搜索区间，避开尖点附近
    all_points <- c(a, zero_points, b)
    for (i in 1:(length(all_points) - 1)) {
      left <- all_points[i] + tol
      right <- all_points[i + 1] - tol
      if (right > left) {
        search_intervals[[length(search_intervals) + 1]] <- c(left, right)
      }
    }
  }

  # 在每个区间寻找导数为零的点
  for (interval in search_intervals) {
    tryCatch({
      # 使用 optimize 寻找极值点
      result <- optimize(fn, interval = interval, maximum = TRUE)

      # 验证是否为极值点（导数接近零）
      if (abs(fn_deriv(result$maximum)) < 0.01) {
        candidate_x <- c(candidate_x, result$maximum)
        candidate_values <- c(candidate_values, result$objective)
      }
    }, error = function(e) {})
  }

  # 去重（保留数值上接近的点中函数值最大的）
  if (length(candidate_x) > 1) {
    unique_indices <- c()
    for (i in 1:length(candidate_x)) {
      if (!any(abs(candidate_x[unique_indices] - candidate_x[i]) < tol)) {
        unique_indices <- c(unique_indices, i)
      }
    }
    candidate_x <- candidate_x[unique_indices]
    candidate_values <- candidate_values[unique_indices]
  }

  # 返回最大值
  max_idx <- which.max(candidate_values)

  return(list(
    maximum = candidate_values[max_idx],
    x_max = candidate_x[max_idx],
    all_candidates = data.frame(x = candidate_x, value = candidate_values)
  ))
}





f_L_SCB <- function(x, cov_theta, mu_star, lambda_best_optim, theta_true){
  x_i <- order_form(x)
  denominator <- sqrt(x_i%*%cov_theta%*%t(t(x_i)))
  lower_bound <- x_i%*%mu_star - lambda_best_optim*denominator
  xTtheta <- x_i%*%t(t(theta_true)) # the true value
  #xTtheta <- x_i%*%t(t(theta_draw)) # the draw theta
  y_l <- xTtheta - lower_bound
  return(y_l)
}

f_U_SCB <- function(x, cov_theta, mu_star, lambda_best_optim, theta_true){
  x_i <- order_form(x)
  denominator <- sqrt(x_i%*%cov_theta%*%t(t(x_i)))
  upper_bound <- x_i%*%mu_star + lambda_best_optim*denominator
  xTtheta <- x_i%*%t(t(theta_true)) # the true value
  #xTtheta <- x_i%*%t(t(theta_draw)) # the draw theta
  y_u <- upper_bound - xTtheta
  return(y_u)
}

#########################################################
###  Huber regression with 5-fold cross-validation  #####
#########################################################
cv_huber <- function(x_scale, Y, k_vals, folds = 5) {
  cv_folds <- caret::createFolds(Y, k = folds, list = TRUE, returnTrain = FALSE) # Create CV folds
  errors <- matrix(NA, nrow = length(k_vals), ncol = folds) # Store results

  for (ii in seq_along(k_vals)) {  # Loop over k values
    k_val <- k_vals[ii]

    for (jj in seq_along(cv_folds)) { # Loop over folds
      test_idx <- cv_folds[[jj]]
      train_idx <- setdiff(seq_along(Y), test_idx)

      X_train <- x_scale[train_idx, , drop = FALSE]  # Split data
      y_train <- Y[train_idx]
      X_test <- x_scale[test_idx, , drop = FALSE]
      y_test <- Y[test_idx]

      huber_fit <- MASS::rlm(y_train ~ ., data = as.data.frame(cbind(y_train, X_train)), psi = MASS::psi.huber, k = k_val) # Fit Huber regression

      y_pred <- predict(huber_fit, newdata = as.data.frame(cbind(rep(1, n=nrow(X_test)),X_test))) # Predict on test data
      errors[ii, jj] <- mean((y_pred - y_test)^2) # Compute mean squared error
    }
  }
  avg_errors <- rowMeans(errors, na.rm = TRUE) # Compute average CV error for each k
  best_k <- k_vals[which.min(avg_errors)] # Get optimal k

  return(list(best_k = best_k, errors = avg_errors))
}



#########################################################
############   Posterior coverage probability  ##########
#########################################################
L_SCB <- function(x, cov_theta, mu_star, lambda_best_optim){
  x_i <- order_form(x)
  denominator <- sqrt(x_i%*%cov_theta%*%t(t(x_i)))
  lower_bound <- x_i%*%mu_star - lambda_best_optim*denominator
  return(lower_bound)
}

U_SCB <- function(x, cov_theta, mu_star, lambda_best_optim){
  x_i <- order_form(x)
  denominator <- sqrt(x_i%*%cov_theta%*%t(t(x_i)))
  upper_bound <- x_i%*%mu_star + lambda_best_optim*denominator
  return(upper_bound)
}

# Function to compute probability at a given x
prob_x_SCB <- function(x, mu_star, cov_theta, dof) {
  x_vec <- order_form(x)
  mu_Y <- sum(x_vec * mu_star)
  sigma_Y <- sqrt(sum(x_vec * (cov_theta %*% x_vec)))
  p <- pt((U_SCB(x) - mu_Y) / sigma_Y, df = dof) - pt((L_SCB(x) - mu_Y) / sigma_Y, df = dof)
  return(p)
}

# Normal density function restricted to [a,b]
truncated_normal_density <- function(x, a, b) {
  denom <- pnorm(b) - pnorm(a)  # Normalization constant
  return(dnorm(x) / denom)
}

# PDF of N(0, 1)
normal_density <- function(x) {
  return(dnorm(x, mean = 0, sd = 1))
}

uniform_density <- function(x, a, b){
  return(dunif(x, min = a, max= b))
}
# Integral function
integrand_SCB <- function(x) {
  if(x_setting == "1"){
    s <- prob_x_SCB(x)*uniform_density(x)
  }else if(x_setting == "2"){
    s <- prob_x_SCB(x)*truncated_normal_density(x)
  }
  return(s)
}

#######    HMC.      #######
# Function to compute probability at a given x for HMC
prob_x_SCB_HMC <- function(x) {

  p <- ecdf_XTtheta(U_SCB(x)) - ecdf_XTtheta(L_SCB(x))
  return(p)
}
# Integral function for HMC
integrand_SCB_HMC <- function(x) {
  if(x_setting == "1"){
    s <- prob_x_SCB_HMC(x)*uniform_density(x)
    #s <- prob_x_SCB_HMC(x)
  }else if(x_setting == "2"){
    s <- prob_x_SCB_HMC(x)*truncated_normal_density(x)
  }
  return(s)
}


find_global_maximum_h_all <- function(a, b, d, cov_mat) {

  # d is a (p+1) vector
  C <- cov_mat                       # (p+1)×(p+1) matrix C is the covariance matrix
  p <- length(d) - 1                 # polynomial order，quadratic: p=2, cubic: p=3

  # polynomial multiplication
  poly_mult <- function(p, q) { # p: order of the polynomial; q: order of another polynomial
    n <- length(p) + length(q) - 1
    result <- numeric(n)
    for (i in seq_along(p))
      result[i:(i + length(q) - 1)] <- result[i:(i + length(q) - 1)] + p[i] * q
    result
  }

  # polynomial derivatives
  poly_deriv <- function(coef) {
    n <- length(coef)
    if (n == 1) return(0)
    coef[-1] * seq_len(n - 1) # coef[-1]: delete the constant term
  }

  # N(x) = (d^T x)^2 = (sum d[i] x^(i-1))^2
  lin_coef <- d
  N_coef <- poly_mult(lin_coef, lin_coef)

  # D(x) = x^T C x = sum_{i,j} C[i,j] x^(i+j-2)
  D_coef <- numeric(2 * p + 1)
  for (i in 1:(p + 1))
    for (j in 1:(p + 1))
      D_coef[i + j - 1] <- D_coef[i + j - 1] + C[i, j]

  # h'(x) = 0  <=>  N'D - ND' = 0
  Np <- poly_deriv(N_coef)
  Dp <- poly_deriv(D_coef)

  NpD <- poly_mult(Np, D_coef)
  NDp <- poly_mult(N_coef, Dp)

  len <- max(length(NpD), length(NDp))
  NpD <- c(NpD, rep(0, len - length(NpD)))
  NDp <- c(NDp, rep(0, len - length(NDp)))
  poly_eq <- NpD - NDp

  # find real roots within (a,b)
  roots <- polyroot(poly_eq)
  real_roots <- Re(roots[abs(Im(roots)) < 1e-6])
  interior_roots <- real_roots[real_roots > a & real_roots < b]


  candidates <- c(a, b, interior_roots)

  # compute T(x) = sqrt(N(x)/D(x))
  T_func <- function(x) {
    xvec <- x^(0:p)
    num  <- as.numeric(d %*% xvec)^2
    den  <- as.numeric(t(xvec) %*% C %*% xvec)
    sqrt(num / den)
  }

  values  <- sapply(candidates, T_func)
  max_idx <- which.max(values)

  return(list(
    maximum        = values[max_idx],
    x_max          = candidates[max_idx],
    all_candidates = data.frame(x = candidates, value = values)
  ))
}

# To compute the critical constant for BSCB and BPCB; To compute PSCP for BSCB and BPCB
sup_T_Bayes_PSCP <- function(a, b, theta_hat, mu_star, cov_mat) {
  d <- theta_hat - as.numeric(mu_star)
  find_global_maximum_h_all(a, b, d = d, cov_mat = cov_mat)
}

## To compute ESCR for BSCB and BPCB
# for BSCB, cov_mat <- cov_theta
# for BPCB, cov_mat <- scale_mat
sup_T_Bayes_ESCR <- function(a, b, theta_true, mu_star, cov_mat) {
  d <- theta_true - as.numeric(mu_star)
  find_global_maximum_h_all(a, b, d = d, cov_mat = cov_mat)
}



## To compute the critical constant for simFSCB
sup_T_simFSCB <- function(a, b, W_sample, cov_mat) { # cov_mat <- XtX_inv
  d <- as.numeric(W_sample)
  find_global_maximum_h_all(a, b, d = d, cov_mat = cov_mat)
}

## To compute the ESCR for FSCB and FPCB
sup_T_Freq_ESCR <- function(a, b, theta_true, lm_theta_hat, cov_mat){ # cov_mat <- (as.numeric(S^2) * inv
  d <- theta_true - as.numeric(lm_theta_hat)
  find_global_maximum_h_all(a, b, d = d, cov_mat = cov_mat)
}





<!-- README.md is generated from README.Rmd. Please edit that file -->

# BSCB

<!-- badges: start -->
<!-- badges: end -->

## Overview

**BSCB** provides methods for constructing Bayesian simultaneous
credible bands (BSCB) and Bayesian pointwise credible bands (BPCB) for
polynomial regression models. The package implements the following
approaches based on different prior specifications:

- **BSCB-C**: Normal-Gamma conjugate prior (empirical Bayes,
  unit-information, or Zellner’s g-prior)
- **BSCB-H**: Non-conjugate prior implemented via HMC
- **BSCB-I-J**: Independent Jeffreys prior (objective Bayesian
  inference)
- **BPCB-I-J**: Bayesian pointwise credible band under the independent
  Jeffreys prior

The methodology is based on the following paper:

- Yang, F., Han, Y., Liu, W., & Hall, I. (2026). [Bayesian simultaneous
  credible bands for polynomial
  regression](https://arxiv.org/abs/2606.28015). arXiv preprint
  arXiv:2606.28015.

A full demo is available
[here](https://github.com/fannyyang73/BSCB/tree/main/demo).

## Installation

``` r
# install.packages("devtools")
devtools::install_github("fannyyang73/BSCB")
```

## Requirements

- R (≥ 4.0.0)
- R packages: `mvtnorm`, `MASS`, `OptimalDesign`, `instantiate`,
  `posterior`

## Quick Start

``` r
library(BSCB)

# Simulate data from a quadratic model
set.seed(123)
n <- 50
x <- seq(-5, 5, length.out = n)
X <- cbind(1, x, x^2)
theta_true <- c(-6, -3, 0.25)
Y <- X %*% theta_true + rnorm(n, sd = 0.2)

# --- BSCB-C: Bayesian simultaneous credible bands under the Normal-Gamma conjugate prior ---
fit_c <- compute_bscb_conjugate(
  X              = X,
  Y              = Y,
  alpha          = 0.05,
  a              = -5,
  b              =  5,
  L              = 500000,
  theta_true     = theta_true,
  hyperparameter = "g_prior",   # "empirical", "unit_info", or "g_prior"
  optimize_type  = "P"          # "P" = polyroot (recommended)
)

# --- BSCB-H: Bayesian simultaneous credible bands under a non-conjugate prior implemented via HMC

mod <- instantiate::stan_package_model(
  name    = "HMC_model",
  package = "BSCB",
  compile = TRUE
)

fit_h <- compute_bscb_hmc(
  X     = X,
  Y     = Y,
  V     = diag(n),
  alpha = alpha,
  a     = a,
  b     = b,
  theta_true = theta_true,
  prior_type = "normal_half_cauchy",
  L     = L,
  draw_num = 10000
)


# --- BSCB-I-J: Bayesian simultaneous credible bands under the Independent Jeffreys prior ---
fit_j <- compute_bscb_ind_jeffreys(
  X     = X,
  Y     = Y,
  alpha = 0.05,
  a     = -5,
  b     =  5,
  L     = 500000,
  theta_true     = theta_true
)

# --- BPCB-I-J: Bayesian pointwise credible bands under the Independent Jeffreys prior ---
fit_p <- compute_bpcb_ind_jeffreys(
  X     = X,
  Y     = Y,
  alpha = 0.05,
  a     = -5,
  b     =  5,
  theta_true     = theta_true
)

# Evaluate bands over a grid and plot
x_seq <- seq(-5, 5, length.out = 500)

plot(x_seq, fit_c$lower_bound(x_seq), type = "l",
     col = "red", lty = 2, lwd = 2,
     ylim = range(c(fit_c$lower_bound(x_seq),
                    fit_c$upper_bound(x_seq), Y)),
     xlab = "x", ylab = "y",
     main = "95% Bayesian Simultaneous Credible Band")
lines(x_seq, fit_c$upper_bound(x_seq), col = "red",  lty = 2, lwd = 2)
lines(x_seq, cbind(1, x_seq, x_seq^2) %*% theta_true,
      col = "blue", lwd = 2)
points(x, Y, pch = 16, col = "gray")
legend("topright",
       legend = c("True curve", "Data", "95% BSCB-C"),
       col    = c("blue", "gray", "red"),
       lty    = c(1, NA, 2),
       pch    = c(NA, 16, NA))
       
# Evaluate coverage
coverage_ESCR(fit_c, optimize_type = "P", verbose = TRUE)
coverage_PSCP(fit_c, draw_num = 10000, optimize_type = "P", verbose = TRUE)
```

## Main Functions

| Function                      | Description                                             |
|-------------------------------|---------------------------------------------------------|
| `compute_bscb_conjugate()`    | BSCB under Normal-Gamma conjugate prior                 |
| `compute_bscb_hmc()`          | BSCB under the non-conjugate prior via HMC              |
| `compute_bscb_ind_jeffreys()` | BSCB under independent Jeffreys prior                   |
| `compute_bpcb_ind_jeffreys()` | BPCB under independent Jeffreys prior                   |
| `coverage_ESCR()`             | Empirical simultaneous coverage rate indicator (0 or 1) |
| `coverage_PSCP()`             | Posterior simultaneous coverage probability estimate    |
| `compute_NG_param()`          | Compute Normal-Gamma posterior parameters               |
| `compute_IJ_param()`          | Compute independent Jeffreys posterior parameters       |

## Key Arguments

| Argument         | Description                           | Options                                                        |
|------------------|---------------------------------------|----------------------------------------------------------------|
| `hyperparameter` | Hyperparameter for Normal-Gamma prior | `"empirical"`, `"unit_info"`, `"g_prior"`                      |
| `optimize_type`  | Method for computing                  | `"P"` (polyroot, recommended), `"G"` (global), `"D"` (DEoptim) |
| `AR_setting`     | Error structure                       | `0` = i.i.d., `1` = AR(1)                                      |
| `L`              | Monte Carlo draws for                 | default `500000`                                               |
| `draw_num`       | Monte Carlo draws for PSCP estimation | default `10000`                                                |

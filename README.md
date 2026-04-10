
<!-- README.md is generated from README.Rmd. Please edit that file -->

# BSCB

<!-- badges: start -->
<!-- badges: end -->

## Overview

**BSCB** provides methods for constructing Bayesian Simultaneous
Credible Bands (BSCB) for polynomial regression models. The package
implements three approaches based on different prior specifications:

- **BSCB-C**: Normal-Gamma conjugate prior (empirical Bayes,
  unit-information, or g-prior)
- **BSCB-J**: Independent Jeffreys prior
- **BSCB-H**: Non-conjugate prior via Hamiltonian Monte Carlo

A full demo is available
[here](https://github.com/fannyyang73/BSCB/tree/main/demo/BSCB_demo.Rmd).

## Installation

``` r
# install.packages("devtools")
devtools::install_github("fannyyang73/BSCB")
```

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

# Fit BSCB with Normal-Gamma conjugate prior
fit <- compute_bscb_conjugate(
  X             = X,
  Y             = Y,
  alpha         = 0.05,       # 1 - 0.05 = 95% credible band
  a             = -5,         # left endpoint of covariate domain
  b             =  5,         # right endpoint of covariate domain
  L             = 50000,      # Monte Carlo draws
  hyperparameter = "g_prior", # "empirical", "unit_info", or "g_prior"
  optimize_type  = "P"        # "P" = polyroot (recommended)
)

# Critical constant
fit$lambda

# Evaluate the band at a grid of x values
x_seq     <- seq(-5, 5, length.out = 500)
lower_vec <- fit$lower_bound(x_seq)
upper_vec <- fit$upper_bound(x_seq)

# Plot
plot(x_seq, lower_vec, type = "l", col = "red", lty = 2, lwd = 2,
     ylim = range(c(lower_vec, upper_vec, Y)),
     xlab = "x", ylab = "y",
     main = "95% Bayesian Simultaneous Credible Band")
lines(x_seq, upper_vec, col = "red", lty = 2, lwd = 2)
lines(x_seq, cbind(1, x_seq, x_seq^2) %*% theta_true, col = "blue", lwd = 2)
points(x, Y, pch = 16, col = "gray")
legend("topright",
       legend = c("True curve", "Data", "95% BSCB"),
       col    = c("blue", "gray", "red"),
       lty    = c(1, NA, 2),
       pch    = c(NA, 16, NA))
```

## Main Functions

| Function                     | Description                                       |
|------------------------------|---------------------------------------------------|
| `compute_bscb_conjugate()`   | BSCB under Normal-Gamma conjugate prior           |
| `compute_bpcb_indJeffreys()` | BSCB under independent Jeffreys prior             |
| `coverage_ESCR()`            | Empirical Simultaneous Coverage Rate (ESCR)       |
| `compute_NG_param()`         | Compute Normal-Gamma posterior parameters         |
| `compute_IJ_param()`         | Compute independent Jeffreys posterior parameters |

## Key Arguments

| Argument         | Description                   | Options                                                        |
|------------------|-------------------------------|----------------------------------------------------------------|
| `hyperparameter` | Hyperparameter specification  | `"empirical"`, `"unit_info"`, `"g_prior"`                      |
| `optimize_type`  | Method for computing sup T(x) | `"P"` (polyroot, recommended), `"G"` (global), `"D"` (DEoptim) |
| `AR_setting`     | Error structure               | `0` = i.i.d., `1` = AR(1)                                      |
| `L`              | Monte Carlo draws for lambda  | default `50000`                                                |

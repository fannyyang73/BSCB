# BSCB 1.0.0

* Initial CRAN submission.

## Features

* Implements six methods for constructing two-sided Bayesian simultaneous 
  credible bands (BSCBs) for the regression curve in univariate polynomial 
  regression over a finite covariate interval:
  - Conjugate Normal-Gamma priors, with empirical Bayes, unit-information, 
    and g-prior hyperparameter specifications
  - Non-conjugate priors fitted via Hamiltonian Monte Carlo (HMC) using 
    'cmdstanr' (Normal-half-Normal and Normal-half-Cauchy priors)
  - A non-informative independent Jeffreys prior approach

* Provides functions for evaluating and comparing method performance:
  - `ESCR()` for computing the empirical simultaneous coverage rate
  - `PSCP()` for computing the posterior simultaneous coverage probability

* Includes a vignette demonstrating usage and comparing methods.

## Documentation

* The methodology is described in Yang, F., Han, Y., Liu, W., & Hall, I. 
  (2026), "Bayesian Simultaneous Credible Bands for Polynomial Regression."

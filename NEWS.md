# BSCB 1.0.2

## Changes

* Added Yang Han, Wei Liu, and Ian Hall as co-authors, reflecting their
  contributions to the methodology described in the accompanying paper
  (arXiv:2606.28015) as well as suggestions on optimize_function.R.


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
  
## Bug fixes

* Fixed an issue where `generate_simulation_data()` with `design_index = 2`
  (D-optimal design) could produce non-reproducible results even when
  `seed` was supplied. The underlying call to `OptimalDesign::od_KL()`
  previously relied on a time-based stopping rule (`t.max`), which made
  the number of completed restarts—and hence the result—dependent on
  machine speed. The search now stops after a fixed number of restarts
  (`rest.max`), ensuring full reproducibility given a fixed `seed`.

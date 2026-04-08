
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Title of the manuscript

<!-- badges: start -->
<!-- badges: end -->

Bayesian Simultaneous Credible Bands for Polynomial Regression

## Authors:

Fei Yang, Yang Han, and Wei Liu

## Configurations:

- R version 4.4.2 (2024-10-31)
- Platform: aarch64-apple-darwin20
- Running under: macOS Ventura 13.7.8
- Running Environment: The Computational Shared Facility (a
  High-Performance Computing cluster at the University of Manchester)

## Execution:

To use the codes for reproducing numerical results, tables and figures
given in the paper and its supplement, please follow the instructions
below.

## Main Manuscript

1.  Illustrative Example - Section 4
    - Reproduce Table 3, Figure 3 and the numerical results in the
      manuscript:
      - Use the file ‘Real_data_example.R’.
2.  Simulation Study – Section 3
    - Reproduce Figure 2 and Table 2 in the manuscript:
      - Use the file ‘Simulation_study.R’.
      - Input the values of confl ($1-\alpha$), Cgamma ($\gamma$), $n$
        and $s$ to compute the values of critical constants $c$ for
        symmetric SCBs and $(c_1,c_2)$ for asymmetric SCBs, and
        $r$-ratios.

``` r
# install.packages("pak")
pak::pak("fannyyang73/BSCB")
```

## Supplementary Materials

3.  Computational Costs
    - Reproduce Tables 1-2 in the supplementary material:
      - Use the file ‘Computational_cost_supple.R’.
4.  Additional Simulation Results
    - Reproduce Tables 3-4 in the supplementary material:
      - Use the file ‘Simulation_study_supple.R’.
      - Input the values of confl ($1-\alpha$), Cgamma ($\gamma$), $n$
        and $s$ to compute the values of critical constants $c$ for
        symmetric SCBs and $(c_1,c_2)$ for asymmetric SCBs, and
        $r$-ratios.

``` r
#library(BSCB)
## basic example code
```

What is special about using `README.Rmd` instead of just `README.md`?
You can include R chunks like so:

``` r
summary(cars)
#>      speed           dist       
#>  Min.   : 4.0   Min.   :  2.00  
#>  1st Qu.:12.0   1st Qu.: 26.00  
#>  Median :15.0   Median : 36.00  
#>  Mean   :15.4   Mean   : 42.98  
#>  3rd Qu.:19.0   3rd Qu.: 56.00  
#>  Max.   :25.0   Max.   :120.00
```

You’ll still need to render `README.Rmd` regularly, to keep `README.md`
up-to-date. `devtools::build_readme()` is handy for this.

You can also embed plots, for example:

<img src="man/figures/README-pressure-1.png" alt="" width="100%" />

In that case, don’t forget to commit and push the resulting figure
files, so they display on GitHub and CRAN.

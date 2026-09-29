
<!-- README.md is generated from README.Rmd. Please edit that file -->

# `citcdf` <img src='man/figures/logo.svg' align="right" height="139" />

[![CRAN
status](https://www.r-pkg.org/badges/version/citcdf)](https://CRAN.R-project.org/package=citcdf)
[![R-CMD-check](https://github.com/sistm/citcdf/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/sistm/citcdf/actions)

## Overview

`citcdf` is a package to perform conditional independence testing using
empirical conditional cumulative distribution function estimations.

The package has two main entry points: `cit_multi()` for gene-wise
conditional independence testing across many outcomes, and `cit_gsa()`
for gene-set analysis. Both accept `test = "asymptotic"` (large samples)
or `test = "permutation"` (small samples). The single-outcome workhorses
`cit_asymp()` and `cit_perm()`, the CCDF estimator `ccdf()`, and the
plotting functions `plot_compare_ccdf()`, `plot.cit_multi()` and
`plot.cit_gsa()` are also exported.

The approach implemented in this package is detailed in the following
article:

> Gauthier M, Agniel D, Thiébaut R & Hejblum BP (2021).
> Distribution-free complex hypothesis testing for single-cell RNA-seq
> differential expression analysis, *bioRxiv*
> [doi:10.1101/2021.05.21.445165](https://doi.org/10.1101/2021.05.21.445165)

## Installation

**`citcdf` is available on
[CRAN](https://CRAN.R-project.org/package=citcdf):**

``` r
install.packages("citcdf")
```

**The development version is available from
[GitHub](https://github.com/sistm/citcdf):**

``` r
# install.packages("remotes")
remotes::install_github("sistm/citcdf")
```

## Example

Here is a basic example which shows how to use `citcdf` with simple
generated data.

``` r
library(citcdf)
```

``` r
## Data Generation
set.seed(123)
n <- 100
X <- data.frame(X1 = as.factor(rbinom(n = n, size = 1, prob = 0.5)))
Y <- replicate(10, (X$X1 == 1) * rnorm(n) + (X$X1 == 0) * rnorm(n, mean = 0.5))
```

``` r
# Hypothesis testing
res_asymp <- cit_multi(M = data.frame(Y = Y), X = X,
                       test = "asymptotic", parallel = FALSE) # asymptotic test
res_asymp$pvals
#>         raw_pval   adj_pval test_statistic
#> Y.1  0.000126502 0.00126502      30.418136
#> Y.2  0.317402043 0.31740204       4.115804
#> Y.3  0.310454844 0.31740204       4.372046
#> Y.4  0.034961503 0.07955292      11.802747
#> Y.5  0.063473425 0.09067632      10.451354
#> Y.6  0.044804258 0.07955292      12.780768
#> Y.7  0.012650228 0.05834245      15.581473
#> Y.8  0.017502734 0.05834245      14.080686
#> Y.9  0.047731750 0.07955292       8.920593
#> Y.10 0.161833735 0.20229217       6.991416
plot(res_asymp)
```

![](man/figures/README-estimation-1.png)<!-- -->

``` r
plot_compare_ccdf(Y=Y[, 1, drop=FALSE], X=X)
```

![](man/figures/README-estimation-2.png)<!-- -->

– Marine Gauthier, Denis Agniel, Sara Fallet, Kalidou Ba, Rodolphe
Thiébaut & Boris Hejblum

*hex illustration by Jérôme Dubois.*

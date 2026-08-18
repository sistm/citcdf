# Estimate the empirical conditional cumulative distribution function

Estimate the empirical conditional cumulative distribution function

## Usage

``` r
ccdf(
  Y,
  X,
  Z = NULL,
  method = c("OLS", "logistic"),
  fast = TRUE,
  space_y = FALSE,
  number_y = 10
)
```

## Arguments

- Y:

  a numeric vector of size `n` containing the preprocessed expressions
  from `n` samples (or cells).

- X:

  a data frame containing numeric or factor vector(s) of size `n`
  containing the variable(s) to be tested (the condition(s) to be
  tested).

- Z:

  a data frame containing numeric or factor vector(s) of size `n`
  containing the covariate(s).

- method:

  a character string indicating which method to use to compute the CCDF,
  either `'OLS'` or `'logistic'`. Default is `'OLS'` for greater
  computational speed.

- fast:

  a logical flag indicating whether the fast implementation of logistic
  regression should be used. Only if `method == 'logistic'`. Default is
  `TRUE`.

- space_y:

  a logical flag indicating whether the y thresholds are spaced. When
  `space_y` is `TRUE`, a regular sequence between the minimum and the
  maximum of the observations is used. Default is `FALSE`.

- number_y:

  an integer value indicating the number of y thresholds (and therefore
  the number of regressions) to perform the test. Only used if `space_y`
  is `TRUE`. Default is `10`.

## Value

A list with the following elements:

- `cdf`: a vector of the cumulative distribution function of a given
  gene.

- `ccdf`: a vector of the conditional cumulative distribution function
  of a given gene, computed given `X`. Only if `Z` is `NULL`.

- `ccdf_nox`: a vector of the conditional cumulative distribution
  function of a given gene, computed given `Z` only (i.e. `X` is
  ignored.). Only if `Z` is not `NULL`.

- `ccdf_x`: a vector of the conditional cumulative distribution function
  of a given gene, computed given `X` and `Z`. Only if `Z` is not
  `NULL`.

- `y_sort`: a vector of the sorted expression points at which the CDF
  and the CCDFs are calculated.

- `x_sort`: a vector of the variables associated with `y_sort`.

- `z_sort`: a vector of the covariates associated with `y_sort`. Only if
  `Z` is not `NULL`.

## Details

The `space_y` / `number_y` grid controls both the resolution of the
statistic and its computational cost. See
[`cit_multi`](https://sistm.github.io/citcdf/reference/cit_multi.md) for
details on this trade-off.

## References

Gauthier M, Agniel D, Thiébaut R & Hejblum BP (2021). Distribution-free
complex hypothesis testing for single-cell RNA-seq differential
expression analysis, *bioRxiv* 445165.
[doi:10.1101/2021.05.21.445165](https://doi.org/10.1101/2021.05.21.445165)
.

## See also

[`plot_compare_ccdf`](https://sistm.github.io/citcdf/reference/plot_compare_ccdf.md),
[`cit_asymp`](https://sistm.github.io/citcdf/reference/cit_asymp.md)

## Examples

``` r

set.seed(123)
n <- 500
X <- as.factor(rbinom(n = n, size = 1, prob = 0.5))
Y <- ((X == 1) * rnorm(n = n, 0, 1)) + ((X == 0) * rnorm(n = n, 0.5, 1))
res <- ccdf(Y, data.frame(X = X), method = "OLS")
```

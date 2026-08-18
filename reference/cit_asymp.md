# Asymptotic test for conditional independence

Test the conditional independence of Y and X given Z.

## Usage

``` r
cit_asymp(Y, X, Z = NULL, space_y = FALSE, number_y = 10, design = NULL)
```

## Arguments

- Y:

  a numeric vector of length `n` to test for conditional independence
  with `X` adjusted on `Z`

- X:

  a data frame of size `n x p` of numeric or factor vector(s) containing
  the variable(s) to be tested for conditional independence against `X`
  adjusted on `Z`.

- Z:

  a data frame of size `n x q` of numeric or factor vector(s) containing
  the covariate(s) to condition the independence test upon.

- space_y:

  a logical flag indicating whether the y thresholds are spaced. When
  `space_y` is `TRUE`, a regular sequence between the minimum and the
  maximum of the observations is used. Default is `FALSE`.

- number_y:

  an integer value indicating the number of y thresholds (and therefore
  the number of regressions) to perform the test. Only used if `space_y`
  is `TRUE`. Default is `10`.

- design:

  an optional (and technical) list of design quantities, as returned by
  the internal `.cit_design(X, Z, n)`. This is used by
  [`cit_multi()`](https://sistm.github.io/citcdf/reference/cit_multi.md),
  to loop-call over many genes while building the model matrix, and
  computing its cross-product and its inverse only once. Default is
  `NULL`, in which case they are computed from `X` and `Z`. Users should
  not be using this argument

## Value

A data frame with the following elements:

- `raw_pval` contains the raw p-values for a given gene.

- `test_statistic` contains the test statistic for a given gene.

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

[`cit_perm`](https://sistm.github.io/citcdf/reference/cit_perm.md),
[`cit_multi`](https://sistm.github.io/citcdf/reference/cit_multi.md),
[`ccdf`](https://sistm.github.io/citcdf/reference/ccdf.md)

## Examples

``` r

set.seed(123)
X <- as.factor(rbinom(n = 100, size = 1, prob = 0.5))
Y <- ((X == 1) * rnorm(n = 100, 0, 1)) + ((X == 0) * rnorm(n = 100, 0.5, 1))
res_asymp <- cit_asymp(Y, data.frame(X = X))


Z <- as.factor(rbinom(n = 100, size = 1, prob = 0.5))
X <- as.numeric(Z) - 1  + rnorm(n = 100, sd = 1)
r <- 500
Y <- replicate(r, as.numeric(Z) - 1)
YY <- (Y == 1) * rnorm(n = 100 * r, 0, 1) + (Y == 0) * rnorm(n = 100 * r, 0.5, 1)
pvals_sim <- sapply(seq_len(r), function(i) {
  cit_asymp(YY[, i], data.frame(X = X), data.frame(Z = Z))$raw_pval
})
hist(pvals_sim) # well calibrated p-values are uniform under the null

quantile(pvals_sim)
#>           0%          25%          50%          75%         100% 
#> 0.0002470792 0.2409864033 0.4958293345 0.7513432132 0.9999476695 
```

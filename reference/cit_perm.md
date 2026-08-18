# Permutation test for conditional independence

Permutation test for conditional independence

## Usage

``` r
cit_perm(
  Y,
  X,
  Z = NULL,
  X_star = NULL,
  n_perm = 100,
  space_y = FALSE,
  number_y = 10
)
```

## Arguments

- Y:

  a numeric vector of length `n` to test for conditional independence
  with `X` adjusted on `Z`

- X:

  a data frame of size `n x p` of numeric or factor vector(s) containing
  the variable(s) to be tested for conditional independence against `X`
  adjusted on `Z`. Multi-variables `X` are supported if `Z` is `NULL`.

- Z:

  a data.frame of size `n x 1` of numeric or factor vector containing
  the covariate to condition the independence test upon. Multiple
  covariates are not supported for permutation.

- X_star:

  a list of `n_perm` permuted designs, as returned by
  [`X_perm`](https://sistm.github.io/citcdf/reference/X_perm.md).
  Default is `NULL`, in which case `X_perm(X, Z, n_perm = n_perm)` is
  called internally. Supply it explicitly when several outcomes must be
  scored against the same permutations, or to avoid redrawing them
  inside a loop; see *Details*.

- n_perm:

  the number of permutations. Default is `100`. When `X_star` is
  supplied it must hold at least `n_perm` elements; only the first
  `n_perm` are used.

- space_y:

  a logical flag indicating whether the y thresholds are spaced out.
  When `space_y` is `TRUE`, a regular sequence between the minimum and
  the maximum of the observations is used. Default is `FALSE`.

- number_y:

  an integer value indicating the number of y thresholds (and therefore
  the number of regressions) to perform the test. Only used if `space_y`
  is `TRUE`. Default is `10`.

## Value

A data frame with the following elements:

- `score` contains the number of permutations whose test statistic is
  greater than or equal to the observed one.

- `raw_pval` contains the raw p-values for a given gene computed from
  `n_perm` permutations.

- `test_statistic` contains the observed test statistic for a given
  gene. It is the same quantity returned by
  [`cit_asymp`](https://sistm.github.io/citcdf/reference/cit_asymp.md).

## Details

The `space_y` / `number_y` grid controls both the resolution of the
statistic and its computational cost. See
[`cit_multi`](https://sistm.github.io/citcdf/reference/cit_multi.md) for
details on this trade-off.

Leaving `X_star` as `NULL` is the convenient form for a single outcome.
Across several outcomes it is not equivalent to supplying one: each call
would draw its own permutations, whereas
[`cit_multi`](https://sistm.github.io/citcdf/reference/cit_multi.md) and
[`cit_gsa`](https://sistm.github.io/citcdf/reference/cit_gsa.md)
deliberately build one pool with
[`X_perm`](https://sistm.github.io/citcdf/reference/X_perm.md) and reuse
it for every gene, so that all genes are scored against the same
permuted designs. Pass a shared `X_star` if you are looping over
outcomes yourself.

## References

Gauthier M, Agniel D, Thiébaut R & Hejblum BP (2021). Distribution-free
complex hypothesis testing for single-cell RNA-seq differential
expression analysis, *bioRxiv* 445165.
[doi:10.1101/2021.05.21.445165](https://doi.org/10.1101/2021.05.21.445165)
.

## See also

[`perm_cont`](https://sistm.github.io/citcdf/reference/perm_cont.md),
[`X_perm`](https://sistm.github.io/citcdf/reference/X_perm.md),
[`cit_multi`](https://sistm.github.io/citcdf/reference/cit_multi.md)

## Examples

``` r

set.seed(123)
X <- data.frame(X = as.factor(rbinom(n = 100, size = 1, prob = 0.5)))
Y <- (X$X == 1) * rnorm(100) + (X$X == 0) * rnorm(100, mean = 0.5)

# the permuted designs are drawn internally when X_star is left NULL
res_perm <- cit_perm(Y, X, n_perm = 10)
res_perm
#>   score   raw_pval test_statistic
#> 1     0 0.09090909        524.969

# supplying them explicitly is equivalent, and is what to do when several
# outcomes must be scored against the same permutations
X_star <- X_perm(X, Z = NULL, n_perm = 10)
res_perm_shared <- cit_perm(Y, X, X_star = X_star, n_perm = 10)

# adjusting for a covariate Z
Z <- data.frame(Z = rnorm(100))
res_perm_adj <- cit_perm(Y, X, Z = Z, n_perm = 10)
res_perm_adj
#>   score   raw_pval test_statistic
#> 1     0 0.09090909        601.138
```

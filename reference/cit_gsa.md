# Conditional independence test for gene set analysis

Conditional independence test for gene set analysis

## Usage

``` r
cit_gsa(
  M,
  X,
  Z = NULL,
  geneset,
  test = c("asymptotic", "permutation"),
  n_perm = 100,
  n_perm_adaptive = c(n_perm, n_perm, n_perm * 3, n_perm * 5),
  thresholds = c(0.1, 0.05, 0.01),
  parallel = interactive(),
  n_cpus = max(1L, detectCores(logical = FALSE) - 1L, na.rm = TRUE),
  adaptive = FALSE,
  space_y = TRUE,
  number_y = 10
)
```

## Arguments

- M:

  a `data.frame` or a `matrix` of size `n x r` containing the different
  Y variables to test for conditional independence with `X` adjusted on
  `Z`.

- X:

  a data frame of size `n x p` of numeric or factor vector(s) containing
  the variable(s) to be tested for conditional independence against `X`
  adjusted on `Z`. Multiple variables (`p>1`) are supported by the
  asymptotic test, and also by the permutation when `Z` is `NULL`.

- Z:

  a data frame of size `n x q` of numeric or factor vector(s) containing
  the covariate(s) to condition the independence test upon. Multiple
  covariates (`q>1`) are only supported by the asymptotic test.

- geneset:

  a vector, a list, a gmt file format or a BiocSet object. If the
  parameter is

  - a vector : corresponds to the gene name of the gene set, must be the
    same as those of the columns of the matrix `M`

  - a list : each elements of the list are a gene set with the names of
    the genes, must be the same as those of the columns of the matrix
    `M`

  - a gmt file format : the genes names of each genes set in the file,
    must be the same as those of the columns of the matrix `M`

  - a BiocSet object : the genes names of each genes set in the object,
    must be the same as those of the columns of the matrix `M`

- test:

  a character string indicating whether the `'asymptotic'` or the
  `'permutation'` test is computed. Default is `'asymptotic'`.

- n_perm:

  the number of permutations. Default is `100`. Only used if
  `test == 'permutation'`.

- n_perm_adaptive:

  a vector of the increasing numbers of adaptive permutations to be
  performed when `adaptive` is `TRUE` if p-values are below
  `thresholds`. `length(n_perm_adaptive)` should be equal to
  `length(thresholds)+1`. Default is
  `c(n_perm, n_perm, n_perm*3, n_perm*5)`.

- thresholds:

  a vector of the decreasing thresholds to compute adaptive permutations
  when `adaptive` is `TRUE`. `length(thresholds)` should be equal to
  `length(n_perm_adaptive)-1`. Default is `c(0.1, 0.05, 0.01)`.

- parallel:

  a logical flag indicating whether parallel computation should be
  enabled. Default is `TRUE` if
  [`interactive()`](https://rdrr.io/r/base/interactive.html) is `TRUE`,
  else is `FALSE`.

- n_cpus:

  an integer indicating the number of cores to be used for the
  computations. Default is
  `max(1L, parallel::detectCores(logical = FALSE) - 1L, na.rm = TRUE)`.
  If `n_cpus = 1`, then sequential computations are used without any
  parallelization.

- adaptive:

  a logical flag indicating whether adaptive permutations should be
  performed. Default is `FALSE` (unlike
  [`cit_multi()`](https://sistm.github.io/citcdf/reference/cit_multi.md)).
  Only used if `test == 'permutation'`.

- space_y:

  a logical flag indicating whether the y thresholds are spaced out.
  When `space_y` is `TRUE`, a regular sequence between the minimum and
  the maximum of the observations is used. If `FALSE`, each unique
  observed expression value is used as a distinct threshold. Default is
  `TRUE`.

- number_y:

  an integer value indicating the number of y thresholds (and therefore
  the number of regressions) to perform the test. Only used if `space_y`
  is `TRUE`. Default is `10`.

## Value

A list with the following elements:

- `which_test`: a character string carrying forward the value of the
  '`test`' argument indicating which test was performed (either
  'asymptotic' or 'permutation').

- `n_perm`: an integer carrying forward the value of the '`n_perm`'
  argument or '`n_perm_adaptive`' indicating the number of permutations
  performed (`NA` if asymptotic test was performed).

- `pvals`: computed p-values. A data frame with one row for each gene
  set, and with 2 columns: the first one '`raw_pval`' contains the raw
  p-values, the second one '`adj_pval`' contains the FDR adjusted
  p-values using Benjamini-Hochberg correction. When
  '`test == "asymptotic"`', a third column '`test_statistic`' contains
  the gene set test statistics. Gene sets with no gene observed in `M`
  yield a warning and `NA` in every column; gene sets only partially
  observed yield a warning and are tested on the measured genes alone.

- `type`: a character string equal to `"gsa"`, identifying the object as
  the result of a gene set analysis.

## Details

The gene-set statistic is the sum of per-gene statistics. For the
permutation test, it is computed with each single permutation of X
shared and applied across all genes in a set (so inter-gene correlation
is preserved).

For the asymptotic test, the null covariance of the stacked threshold
indicators is estimated empirically (`crossprod(temp) / n`). The closed
form \\min(p_i, p_j) - p_i p_j\\ used by
[`cit_asymp`](https://sistm.github.io/citcdf/reference/cit_asymp.md)
only holds within a gene (because the product of two threshold
indicators of the same `Y` is itself an indicator). Across two genes,
that same expectation is their joint distribution function, which the
marginal proportions do not determine as there is inter-gene correlation
present. The empirical estimator coincides with the closed form on the
within-gene diagonal blocks, and additionally supplies the between-gene
blocks which carry the inter-gene correlation needed by the summed
gene-set statistic.

The `space_y` / `number_y` grid controls both the resolution of the
statistic and its computational cost. See
[`cit_multi`](https://sistm.github.io/citcdf/reference/cit_multi.md) for
details on this trade-off.

## See also

[`cit_multi`](https://sistm.github.io/citcdf/reference/cit_multi.md),
[`plot.cit_gsa`](https://sistm.github.io/citcdf/reference/plot.cit_gsa.md)

## Examples

``` r
# Two conditions and 30 genes, split into two sets: in "responder" every of
# the 10 genes shifts slightly with X, in "null" none of the remaining 20 does.
set.seed(123)
n <- 100
X <- data.frame(X = as.factor(rbinom(n, size = 1, prob = 0.5)))
M <- matrix(rnorm(n * 30), nrow = n, dimnames = list(NULL, paste0("g", 1:30)))
M[, 1:10] <- M[, 1:10] + 0.3 * (as.numeric(X$X) - 1)
geneset <- list(responder = paste0("g", 1:10), null = paste0("g", 11:30))

res <- cit_gsa(M = M, X = X, geneset = geneset,
  test = "asymptotic", parallel = FALSE)
res$pvals
#>            raw_pval   adj_pval test_statistic
#> responder 0.0199274 0.03985479       71.07193
#> null      0.8877189 0.88771885       60.42345

# Single gene shifts are too small to be detected on their own,
# but together the set is significant.
per_gene <- cit_multi(M = as.data.frame(M[, 1:10]), X = X,
  test = "asymptotic", parallel = FALSE)
min(per_gene$pvals$adj_pval)  # no single gene survives the correction
#> [1] 0.1505369
res$pvals["responder", ]      # the set does
#>            raw_pval   adj_pval test_statistic
#> responder 0.0199274 0.03985479       71.07193

# a single gene set may be given as a plain character vector of M colnames
cit_gsa(M = M, X = X, geneset = paste0("g", 1:10),
  test = "asymptotic", parallel = FALSE)$pvals
#>    raw_pval  adj_pval test_statistic
#> 1 0.0199274 0.0199274       71.07193

# adjusting for a covariate
Z <- data.frame(Z = rnorm(n))
cit_gsa(M = M, X = X, Z = Z, geneset = geneset,
  test = "asymptotic", parallel = FALSE)$pvals
#>             raw_pval   adj_pval test_statistic
#> responder 0.01583112 0.03166225       73.69177
#> null      0.88263643 0.88263643       61.42700

# \donttest{
# The permutation test applies each single permutation of X to every gene of a
# set at once, so the correlation between genes is carried into the null.
cit_gsa(M = M, X = X, geneset = geneset,
  test = "permutation", n_perm = 100, parallel = FALSE)$pvals
#> Computing 100 permutations...
#>             raw_pval   adj_pval test_statistic
#> responder 0.01980198 0.03960396       71.07193
#> null      0.96039604 0.96039604       60.42345

# genes listed in a set but absent from M are dropped, with a warning
cit_gsa(M = M, X = X,
  geneset = list(partly_measured = c(paste0("g", 1:5), "absent1")),
  test = "asymptotic", parallel = FALSE)$pvals
#> Warning:  Some genes from geneset partly_measured are not observed in expression data
#>                  raw_pval  adj_pval test_statistic
#> partly_measured 0.1075917 0.1075917       30.46563
# }
```

# Permutation procedure when Z is continuous

Permutation procedure when Z is continuous

## Usage

``` r
perm_cont(X, Z)
```

## Arguments

- X:

  a numeric or factor vector of length `n` containing the variable to be
  tested (the condition to be tested).

- Z:

  a numeric vector of length `n` containing the covariate. Multiple
  variables are not allowed.

## Value

`X_star` a vector of permuted `X`.

## Details

`X_star` is a conditional permutation of `X`: draws are without
replacement, so `X_star` is a bijective rearrangement of `X` and its
marginal distribution is preserved exactly. The permut for observation
`i` is drawn with a weight that decreases with the distance between
fitted values of `X` given `Z`, so they are matched on `Z`.

For binary `X` the fitted value is a linear-probability approximation of
the propensity score (matching on the propensity score balances `Z`,
\[Rosenbaum and Rubin 1983\] and allows to match on that single number
rather than on `Z`). Of note, this approximation can fall outside
`[0, 1]` when the association between `X` and `Z` is strong.

This is the conditional permutation test (CPT) of Berrett et al. (2020),
which draws a permutation with probability proportional to the
likelihood of the permuted assignment under the conditional law of `X`
given `Z`. The sequential draw used here approximates that distribution
rather than sampling from it exactly. Permuting rather than resampling
`X`, as the conditional randomisation test (CRT) of Candes et al. (2018)
does, preserves the empirical marginal distribution of `X` exactly
regardless of the error in the fitted conditional law model (and never
produces degenerate single-level designs).

Permutation weights are Gaussian in the distance between fitted values,
with bandwidth `sd(fit) * n^(-1/3)` (polynomial weights leave too much
mass on distant candidates, so the neighbourhood was not local and the
`X`-`Z` relationship was not preserved). Of note, the bandwith exponent
is larger than Silverman's `1/5` because his rate is optimized for
density estimation but leaves a first-order bias that invalidates
inference built on it (Hall, 1992; Armstrong and Kolesar, 2020);
`o(n^(-1/4))` is the condition Kim et al. (2022) prove for the analogous
local permutation test.

## References

Berrett TB, Wang Y, Barber RF, Samworth RJ (2020). The conditional
permutation test for independence while controlling for confounders.
*Journal of the Royal Statistical Society Series B*, **82**(1), 175-197.
[doi:10.1111/rssb.12340](https://doi.org/10.1111/rssb.12340)

Candes E, Fan Y, Janson L, Lv J (2018). Panning for gold: 'model-X'
knockoffs for high dimensional controlled variable selection. *Journal
of the Royal Statistical Society Series B*, **80**(3), 551-577.
[doi:10.1111/rssb.12265](https://doi.org/10.1111/rssb.12265)

Hemerik J, Goeman JJ (2018). Exact testing with random permutations.
*TEST*, **27**(4), 811-825.
[doi:10.1007/s11749-017-0571-1](https://doi.org/10.1007/s11749-017-0571-1)

Rosenbaum PR, Rubin DB (1983). The central role of the propensity score
in observational studies for causal effects. *Biometrika*, **70**(1),
41-55.
[doi:10.1093/biomet/70.1.41](https://doi.org/10.1093/biomet/70.1.41)

Kim I, Neykov M, Balakrishnan S, Wasserman L (2022). Local permutation
tests for conditional independence. *The Annals of Statistics*,
**50**(6), 3388-3414.
[doi:10.1214/22-AOS2233](https://doi.org/10.1214/22-AOS2233)

Hall P (1992). Effect of bias estimation on coverage accuracy of
bootstrap confidence intervals for a probability density. *The Annals of
Statistics*, **20**(2), 675-694.
[doi:10.1214/aos/1176348651](https://doi.org/10.1214/aos/1176348651)

Armstrong TB, Kolesar M (2020). Simple and honest confidence intervals
in nonparametric regression. *Quantitative Economics*, **11**(1), 1-39.
[doi:10.3982/QE1199](https://doi.org/10.3982/QE1199)

## See also

[`X_perm`](https://sistm.github.io/citcdf/reference/X_perm.md)

## Examples

``` r

set.seed(123)
X <- rbinom(n = 100, size = 1, prob = 0.5)
Z <- rnorm(100, 0, 1)
X_star <- perm_cont(X, Z)
table(X, X_star)
#>    X_star
#> X    0  1
#>   0 30 23
#>   1 23 24
```

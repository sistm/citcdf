# citcdf: Conditional Independence Testing with Cumulative Distribution Functions

Distribution-free conditional independence testing built on estimates of
the conditional cumulative distribution function (CCDF).

## Main functions

- [`cit_multi`](https://sistm.github.io/citcdf/reference/cit_multi.md):
  gene-wise testing across many outcomes.

- [`cit_gsa`](https://sistm.github.io/citcdf/reference/cit_gsa.md):
  gene-set analysis.

- [`cit_asymp`](https://sistm.github.io/citcdf/reference/cit_asymp.md) /
  [`cit_perm`](https://sistm.github.io/citcdf/reference/cit_perm.md):
  the single-outcome asymptotic and permutation tests.

- [`ccdf`](https://sistm.github.io/citcdf/reference/ccdf.md): the CCDF
  estimator the tests are built on.

- [`plot_compare_ccdf`](https://sistm.github.io/citcdf/reference/plot_compare_ccdf.md):
  diagnostic CCDF plots.

## Note on defaults

[`cit_multi()`](https://sistm.github.io/citcdf/reference/cit_multi.md)
and [`cit_gsa()`](https://sistm.github.io/citcdf/reference/cit_gsa.md)
use `space_y = TRUE` with `number_y = 10` (for computational speed),
while [`ccdf()`](https://sistm.github.io/citcdf/reference/ccdf.md),
[`cit_asymp()`](https://sistm.github.io/citcdf/reference/cit_asymp.md)
and [`cit_perm()`](https://sistm.github.io/citcdf/reference/cit_perm.md)
default to `space_y = FALSE`, i.e. every distinct observed value is a
threshold.

## References

Gauthier M, Agniel D, Thiébaut R & Hejblum BP (2021). Distribution-free
complex hypothesis testing for single-cell RNA-seq differential
expression analysis, *bioRxiv* 445165.
[doi:10.1101/2021.05.21.445165](https://doi.org/10.1101/2021.05.21.445165)

## See also

Useful links:

- <https://github.com/sistm/citcdf>

- <https://sistm.github.io/citcdf/>

- Report bugs at <https://github.com/sistm/citcdf/issues>

## Author

**Maintainer**: Boris P. Hejblum <boris.hejblum@u-bordeaux.fr>

Authors:

- Boris P. Hejblum <boris.hejblum@u-bordeaux.fr>

- Denis Agniel <denis.agniel@gmail.com>

- Sara Fallet <sara.fallet@u-bordeaux.fr>

- Marine Gauthier <marine.gauthier@epoch-intelligence.fr>

Other contributors:

- Kalidou Ba <kalidou.ba@u-bordeaux.fr> \[contributor\]

- Pierre Neuvial <pierre.neuvial@math.univ-toulouse.fr> \[contributor\]

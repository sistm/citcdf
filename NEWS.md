# citcdf 1.1.0
 * `space_y = TRUE` now uses the same y-threshold grid across all tests: `cit_asymp()` and `cit_gsa()` adopt `cit_perm()`'s upper endpoint `max(Y)` instead of the second-largest observation. Asymptotic results will be different than previously with defaults, although calibration remains unchanged.
 * first CRAN release

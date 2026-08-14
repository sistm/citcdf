# citcdf 1.1.0
 
 * `space_y = TRUE` now uses the same y-threshold grid across all tests: 
 `cit_asymp()`, and `cit_gsa()` adopt `cit_perm()`'s upper endpoint `max(Y)` 
 instead of the second-largest observation. Asymptotic results will be different 
 than previously with defaults, although calibration remains unchanged.
 
  * similarly, `ccdf()` now uses the same y-threshold grid as the tests when 
  `space_y = TRUE` (starting at the smallest non-zero observation rather than 
  the second smallest). CCDF values from `ccdf(space_y = TRUE)` therefore differ 
  from previous versions (the default `space_y = FALSE` is unaffected). It also 
  no longer errors when all `Y` are equal.
 
 * first CRAN release

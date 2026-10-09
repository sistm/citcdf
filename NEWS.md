# citcdf 1.2.1

* **Covariance** in asymptotic tests (`cit_asymp()`, `cit_multi()`, `cit_gsa()`) 
  is now estimated with a *heteroskedasticity-robust sandwich estimator* by default. 
  Former factorized estimator was only valid if `Y` is independent of both
  `X` and `Z`, and was too conservative when `Z` actually affects `Y`. Test 
  statistics are unchanged but asymptotic p-values change. A new argument `variance`   
  is either `"sandwich"` (the new default) or `"independent"` (older
  estimator for backward compatibility with citcdf 1.1.0 results).

* New argument `residuals`: `"full"` (default, like before) or `"restricted"` 
  (null-model, score-type, more conservative) residuals in the sandwich estimator.
  
* New argument `small_sample_corr` in asymptotic tests, correcting the
  *anti-conservativeness* of the sandwich estimator in small samples or with
  many covariates: `"HC2"` (residuals divided by `sqrt(1 - leverage)`), `"HC1"`
  (estimator multiplied by `n/(n-d)`) or `"none"` (as in 1.2.0). **Default is
  `"HC2"` in `cit_asymp()` and `cit_multi()`** and `"none"` in `cit_gsa()`.



# citcdf 1.1.0
 
 * `perm_cont()`:
     - now draws **without replacement**, so the conditional permutation 
     `X_star` is a genuine permutation of `X` (rather than a re-sample) and 
     preserves `X`'s empirical marginal distribution. 
     Because this changes both the RNG stream and the null distribution, permutation 
     p-values from `cit_perm()`,`cit_multi()` and `cit_gsa()` will differ from 
     previous versions. Of note, previous scheme could be anti-conservative, and 
     some limiting cases were badly handled.
     - now uses QR decomposition through `.lm.fit()` to be robust to badly scaled 
     covariates (eg library sizes).
     - now weights draws with a **Gaussian kernel** on the fitted
     values (bandwidth `sd(fit) * n^(-1/3)`) rather than previous
     `1/(difference in fitted values)^2` (each observation effectively drew 
     from its two nearest neighbours whatever the sample size, leading to bad 
     calibration of the test with too many false positives). This also changes
     the null distribution.
 
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

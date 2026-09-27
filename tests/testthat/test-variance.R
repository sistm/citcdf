# variance = "independent": deprecated factorised estimator, kept for
# backward compatibility. It must reproduce citcdf 1.1.0 exactly.

snap_fixture_110 <- function() {                  # = snap_fixture() in test-snapshot.R
  set.seed(20240301)
  n <- 60
  X <- data.frame(X = as.factor(rbinom(n, 1, 0.5)))
  Z <- data.frame(Z = rnorm(n))
  Y <- 0.5 * Z$Z + 1.2 * (X$X == 1) + rnorm(n)
  M <- cbind(Y, replicate(4, 0.5 * Z$Z + rnorm(n)))
  colnames(M) <- paste0("g", seq_len(ncol(M)))
  list(X = X, Z = Z, Y = Y, M = M)
}

test_that("variance = 'independent' reproduces the citcdf 1.1.0 p-values", {
  f <- snap_fixture_110()
  indep <- function(expr) suppressWarnings(expr)
  # reference values: tests/testthat/_snaps/snapshot.md at tag v1.1.0
  # (relative tolerance 1e-6: the saddlepoint approximation is sensitive to the
  # BLAS in its last digits when Sigma has numerically null eigenvalues)
  expect_equal(indep(cit_asymp(f$Y, f$X, variance = "independent"))$raw_pval,
    2.96693654018144e-05, tolerance = 1e-6)
  expect_equal(indep(cit_asymp(f$Y, f$X, f$Z, variance = "independent"))$raw_pval,
    4.18507499990775e-05, tolerance = 1e-6)
  expect_equal(indep(cit_asymp(f$Y, f$X, space_y = TRUE, number_y = 10,
    variance = "independent"))$raw_pval, 5.029411060088e-05, tolerance = 1e-6)
  multi <- indep(cit_multi(M = data.frame(f$M), X = f$X, Z = f$Z,
    test = "asymptotic", parallel = FALSE, variance = "independent"))$pvals
  expect_equal(multi$raw_pval, c(7.2869577061517e-05, 0.409146461282488,
    0.490150214439008, 0.610306105045646, 0.865286501017491), tolerance = 1e-6)
  gsa <- indep(cit_gsa(M = data.frame(f$M), X = f$X,
    geneset = list(a = c("g1", "g2"), b = c("g3", "g4", "g5")),
    test = "asymptotic", parallel = FALSE, variance = "independent"))$pvals
  expect_equal(gsa$raw_pval, c(0.000689771136486358, 0.859197717333256),
    tolerance = 1e-6)
})

test_that("variance = 'independent' warns once per call, the default does not", {
  f <- snap_fixture_110()
  n_depr <- function(expr) sum(grepl("deprecated", testthat::capture_warnings(expr)))
  expect_equal(n_depr(cit_asymp(f$Y, f$X, variance = "independent")), 1)
  expect_equal(n_depr(cit_multi(M = data.frame(f$M), X = f$X, test = "asymptotic",
    parallel = FALSE, variance = "independent")), 1)
  expect_equal(n_depr(cit_gsa(M = data.frame(f$M), X = f$X,
    geneset = list(a = c("g1", "g2"), b = c("g3", "g4")),
    parallel = FALSE, variance = "independent")), 1)
  expect_no_warning(cit_asymp(f$Y, f$X, f$Z))
  expect_error(cit_asymp(f$Y, f$X, variance = "foo"))
})

test_that("without Z, both estimators agree asymptotically under the null", {
  set.seed(8)
  n <- 3000
  X <- data.frame(X = as.factor(rbinom(n, 1, 0.5)))
  Y <- rnorm(n)
  d <- .cit_design(X, NULL, n)
  y <- .cit_y_grid(Y, TRUE, 10)
  D <- outer(Y, y[-length(y)], "<=") * 1
  expect_equal(sum(.cit_sandwich_ev(D, d)), sum(.cit_independent_ev(D, d)),
    tolerance = 0.05)
})

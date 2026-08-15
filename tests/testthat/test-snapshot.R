# Numerical snapshots. 
# 
# These pin the VALUES, not just the shape, so that
# a refactor meant to be numerically neutral can prove it.
#
# The references live in _snaps/snapshot.json and are compared with a tolerance,
# not bit-for-bit: several of those changes are expected to move results by
# ~1e-14, and cancelling the self-cancelling outer products in cit_gsa() is not
# bit-exact even in principle.
#
# If a numerical change is DELIBERATE, run testthat::snapshot_accept("snapshot")
# and say why in NEWS.md. If it is not, this file is what tells you.
#
# expect_snapshot_value() skips on CRAN by default: the last
# digits depend on the BLAS, so pinning them is useful against our own refactors
# but a poor reason to fail somebody else's check machine.

snap_fixture <- function() {
  set.seed(20240301)
  n <- 60
  X <- data.frame(X = as.factor(rbinom(n, 1, 0.5)))
  Z <- data.frame(Z = rnorm(n))
  Y <- 0.5 * Z$Z + 1.2 * (X$X == 1) + rnorm(n)
  M <- cbind(Y, replicate(4, 0.5 * Z$Z + rnorm(n)))
  colnames(M) <- paste0("g", seq_len(ncol(M)))
  list(X = X, Z = Z, Y = Y, M = M)
}

test_that("ccdf estimates are numerically stable", {
  f <- snap_fixture()
  expect_snapshot_value(tolerance = 1e-10, style = "deparse",
    list(ols      = as.numeric(ccdf(Y = f$Y, X = f$X)$ccdf),
         ols_z    = as.numeric(ccdf(Y = f$Y, X = f$X, Z = f$Z)$ccdf_x),
         logistic = as.numeric(ccdf(Y = f$Y, X = f$X, method = "logistic")$ccdf)))
})

test_that("test statistics and p-values are numerically stable", {
  f <- snap_fixture()
  asymp <- function(...) unlist(cit_asymp(...)[c("test_statistic", "raw_pval")])
  multi <- cit_multi(M = data.frame(f$M), X = f$X, Z = f$Z,
    test = "asymptotic", parallel = FALSE)$pvals
  gsa <- cit_gsa(M = data.frame(f$M), X = f$X,
    geneset = list(a = c("g1", "g2"), b = c("g3", "g4", "g5")),
    test = "asymptotic", parallel = FALSE)$pvals
  expect_snapshot_value(tolerance = 1e-10, style = "deparse",
    list(asymp      = as.numeric(asymp(Y = f$Y, X = f$X)),
         asymp_z    = as.numeric(asymp(Y = f$Y, X = f$X, Z = f$Z)),
         asymp_grid = as.numeric(asymp(Y = f$Y, X = f$X, space_y = TRUE,
                                       number_y = 10)),
         multi      = as.numeric(c(multi$test_statistic, multi$raw_pval)),
         gsa        = as.numeric(c(gsa$test_statistic, gsa$raw_pval))))
})

test_that("results are reproducible", {
  f <- snap_fixture()
  asymp <- function() cit_multi(M = data.frame(f$M), X = f$X, Z = f$Z,
    test = "asymptotic", parallel = FALSE)$pvals$test_statistic
  expect_identical(asymp(), asymp())
  perm <- function() {
    set.seed(99)
    cit_multi(M = data.frame(f$M), X = f$X, Z = f$Z, test = "permutation",
      n_perm = 30, adaptive = FALSE, parallel = FALSE)$pvals$raw_pval
  }
  expect_identical(perm(), perm())
})

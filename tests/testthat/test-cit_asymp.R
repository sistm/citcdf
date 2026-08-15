test_that("cit_asymp returns a valid test across designs, and has power", {
  set.seed(1)
  n <- 200
  Z <- data.frame(Z = rnorm(n))
  designs <- list(
    binary   = data.frame(X = as.factor(rbinom(n, 1, 0.5))),
    twocol   = data.frame(X1 = as.factor(rbinom(n, 1, 0.5)),   # ncol(X) >= 2:
      X2 = as.factor(rbinom(n, 1, 0.5))),  # Sigma was non-PSD here
    threelev = data.frame(X = as.factor(sample(0:2, n, TRUE))) # nrow(H) > 1
  )
  for (nm in names(designs)) {
    X <- designs[[nm]]
    for (zz in list(NULL, Z)) {
      res <- cit_asymp(Y = 0.8 * Z$Z + rnorm(n), X = X, Z = zz)
      expect_true(all(c("raw_pval", "test_statistic") %in% names(res)), info = nm)
      expect_true(is.finite(res$test_statistic), info = nm)
      expect_true(res$raw_pval >= 0 && res$raw_pval <= 1, info = nm)
    }
  }
  X <- designs$binary
  expect_lt(cit_asymp(Y = rnorm(n) + 1.5 * (X$X == 1), X = X)$raw_pval, 0.01)
  expect_gt(cit_asymp(Y = rnorm(n), X = X)$raw_pval, 0.05)
})

test_that("a finer space_y grid converges to the full grid, and the row order is irrelevant", {
  set.seed(2)
  n <- 120
  X <- data.frame(X = as.factor(rbinom(n, 1, 0.5)))
  Y <- rnorm(n) + 1.2 * (X$X == 1)
  full <- cit_asymp(Y = Y, X = X, space_y = FALSE)
  coarse <- cit_asymp(Y = Y, X = X, space_y = TRUE, number_y = 5)$test_statistic
  finer <- cit_asymp(Y = Y, X = X, space_y = TRUE, number_y = 50)$test_statistic
  expect_lt(abs(finer - full$test_statistic), abs(coarse - full$test_statistic))
  o <- sample.int(n)
  expect_equal(cit_asymp(Y = Y[o], X = X[o, , drop = FALSE])$test_statistic,
    full$test_statistic)
})

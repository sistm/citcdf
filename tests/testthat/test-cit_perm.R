test_that("X_perm returns a list of permuted designs of the requested length", {
  set.seed(2)
  X <- data.frame(X = as.factor(rbinom(50, 1, 0.5)))
  xs <- X_perm(X, Z = NULL, n_perm = 12)
  expect_type(xs, "list")
  expect_length(xs, 12)
  expect_s3_class(xs[[1]], "data.frame")
  expect_equal(nrow(xs[[1]]), nrow(X))
})

test_that("cit_perm works with an X_star pool", {
  set.seed(3)
  X <- data.frame(X = as.factor(rbinom(80, 1, 0.5)))
  Y <- (X$X == 1) * rnorm(80) + (X$X == 0) * rnorm(80, mean = 0.5)
  xs <- X_perm(X, Z = NULL, n_perm = 20)
  res <- cit_perm(Y, X, X_star = xs, n_perm = 20)
  expect_true(res$raw_pval >= 0 && res$raw_pval <= 1)
})

test_that("cit_perm returns the observed statistic alongside score and p-value", {
  set.seed(5)
  n <- 100
  X <- data.frame(X = as.factor(rbinom(n, 1, .5)))
  Y <- rnorm(n) + 0.6 * (as.numeric(X$X) - 1)
  Xs <- X_perm(X, NULL, n_perm = 30)
  res <- cit_perm(Y, X, NULL, X_star = Xs, n_perm = 30, space_y = FALSE)

  expect_named(res, c("score", "raw_pval", "test_statistic"))
  expect_true(is.finite(res$test_statistic) && res$test_statistic >= 0)
  # score is a count bounded by n_perm -- NOT the statistic
  expect_true(res$score >= 0 && res$score <= 30)
})

test_that("cit_asymp and cit_perm report the SAME observed statistic", {
  set.seed(5); n <- 100
  X <- data.frame(X = as.factor(rbinom(n, 1, .5)))
  Y <- rnorm(n) + 0.6 * (as.numeric(X$X) - 1)
  Xs <- X_perm(X, NULL, n_perm = 30)
  
  # space_y = FALSE so both use sort(unique(Y)); see C5 before relaxing this
  a <- cit_asymp(Y, X, space_y = FALSE)
  p <- cit_perm(Y, X, NULL, X_star = Xs, n_perm = 30, space_y = FALSE)
  
  expect_identical(a$test_statistic, p$test_statistic)   # verified: diff exactly 0
})

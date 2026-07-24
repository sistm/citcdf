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
  expect_true("raw_pval" %in% colnames(res))
  expect_true(res$raw_pval >= 0 && res$raw_pval <= 1)
})

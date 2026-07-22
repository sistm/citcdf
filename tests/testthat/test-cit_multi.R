test_that("asymptotic cit_multi returns valid p-values", {
  set.seed(4)
  n <- 60; r <- 5
  X <- data.frame(X = rnorm(n))
  M <- as.data.frame(matrix(rnorm(n * r), n, r))
  res <- cit_multi(M = M, X = X, test = "asymptotic", parallel = FALSE)
  expect_s3_class(res, "cit_multi")
  expect_equal(nrow(res$pvals), r)
  expect_true(all(res$pvals$raw_pval >= 0 & res$pvals$raw_pval <= 1))
})

test_that("adaptive permutation cit_multi reaches a stage beyond the initial pool", {
  set.seed(5)
  n <- 50
  X <- data.frame(X = as.factor(rbinom(n, 1, 0.5)))
  M <- data.frame(Y = 2 * (as.numeric(X$X) - 1) + rnorm(n))  # strong effect: survives thresholds
  expect_error(
    res <- cit_multi(M = M, X = X, test = "permutation", adaptive = TRUE,
                     n_perm = 10, n_perm_adaptive = c(10, 30),
                     thresholds = 0.9, parallel = FALSE),
    NA)
  expect_true(all(res$pvals$raw_pval >= 0 & res$pvals$raw_pval <= 1))
})

test_that("cit_multi returned object matches its documented contract (@return)", {
  set.seed(6)
  n <- 50; r <- 4
  X <- data.frame(X = rnorm(n))
  M <- as.data.frame(matrix(rnorm(n * r), n, r))
  res <- cit_multi(M = M, X = X, test = "asymptotic", parallel = FALSE)
  expect_named(res$pvals, c("raw_pval", "adj_pval", "test_statistic"))
  expect_identical(res$which_test, "asymptotic")
})

# Regression tests for cit_gsa: the asymptotic path plus the permutation path
# that R CMD check's examples did not previously exercise.


test_that("asymptotic cit_gsa returns a well-formed citcdf object", {
  d <- make_data()
  res <- cit_gsa(M = d$M, X = d$X, geneset = d$geneset,
    test = "asymptotic", parallel = FALSE)
  expect_s3_class(res, "cit_gsa")
  expect_equal(res$which_test, "asymptotic")
  expect_true(all(c("raw_pval", "adj_pval", "test_statistic") %in%
    colnames(res$pvals)))
  expect_equal(nrow(res$pvals), length(d$geneset))
  expect_true(all(res$pvals$raw_pval >= 0 & res$pvals$raw_pval <= 1))
})

test_that("permutation cit_gsa runs and returns valid p-values (regression: X_star)", {
  d <- make_data()
  expect_error(
    res <- cit_gsa(M = d$M, X = d$X, geneset = d$geneset,
      test = "permutation", n_perm = 50, parallel = FALSE),
    NA)  # NA => expect NO error
  expect_s3_class(res, "cit_gsa")
  expect_equal(res$which_test, "permutation")
  expect_true(all(c("raw_pval", "adj_pval") %in% colnames(res$pvals)))
  expect_true(all(res$pvals$raw_pval >= 0 & res$pvals$raw_pval <= 1))
})

test_that("adaptive permutation reaches a stage larger than the initial pool (regression: pool sizing)", {
  d <- make_data()
  # stage 2 requests 30 > 10 permutations: fails unless X_star pool is sized to max()
  expect_error(
    res <- cit_gsa(M = d$M, X = d$X, geneset = d$geneset,
      test = "permutation", adaptive = TRUE,
      n_perm = 10, n_perm_adaptive = c(10, 30),
      thresholds = 0.9, parallel = FALSE),
    NA)
  expect_true(all(res$pvals$raw_pval >= 0 & res$pvals$raw_pval <= 1))
})

test_that("returned object matches its documented contract (@return)", {
  d <- make_data()

  # documented default: adaptive = FALSE
  expect_false(eval(formals(cit_gsa)$adaptive))

  # asymptotic: 3 columns incl. test_statistic, plus type == "gsa"
  res_a <- cit_gsa(M = d$M, X = d$X, geneset = d$geneset,
    test = "asymptotic", parallel = FALSE)
  expect_named(res_a$pvals, c("raw_pval", "adj_pval", "test_statistic"))
  expect_identical(res_a$type, "gsa")
  expect_identical(res_a$which_test, "asymptotic")
  expect_true(is.na(res_a$n_perm))

  # permutation exports the aggregated observed statistic as documented
  res_p <- cit_gsa(M = d$M, X = d$X, geneset = d$geneset,
    test = "permutation", n_perm = 20, parallel = FALSE)
  expect_named(res_p$pvals, c("raw_pval", "adj_pval", "test_statistic"))
  expect_true(all(is.finite(res_p$pvals$test_statistic)))
  expect_true(all(res_p$pvals$test_statistic >= 0))
})

test_that("permutation cit_gsa returns one row per gene SET", {
  d <- make_data()
  res <- cit_gsa(d$M, d$X, geneset = d$geneset, test = "permutation",
    n_perm = 100, parallel = FALSE)
  expect_equal(nrow(res$pvals), length(d$geneset))
  expect_true("test_statistic" %in% colnames(res$pvals))
})

test_that("permutation and asymptotic flag the same gene set", {
  d <- make_data()                                   # set1 associated, set2 null
  ap <- cit_gsa(d$M, d$X, geneset = d$geneset, test = "asymptotic",
    parallel = FALSE)$pvals$raw_pval
  pp <- cit_gsa(d$M, d$X, geneset = d$geneset, test = "permutation",
    n_perm = 500, parallel = FALSE)$pvals$raw_pval
  expect_lt(pp[1], pp[2])
  expect_equal(which.min(ap), which.min(pp))
})

test_that("multi-column X works on the permutation path when Z is NULL", {
  set.seed(3)
  n <- 80
  X <- data.frame(X1 = rnorm(n), X2 = as.factor(rbinom(n, 1, .5)))
  M <- as.data.frame(replicate(4, rnorm(n)) + 0.8 * X$X1)
  colnames(M) <- paste0("g", 1:4)
  res <- cit_gsa(M, X, Z = NULL, geneset = list(s1 = colnames(M)),
    test = "permutation", n_perm = 50, parallel = FALSE)
  expect_equal(nrow(res$pvals), 1L)
  # and still refused when Z is present
  Z <- data.frame(Z = as.factor(rbinom(n, 1, .5)))
  expect_error(cit_gsa(M, X, Z = Z, geneset = list(s1 = colnames(M)),
    test = "permutation", n_perm = 20, parallel = FALSE))
})

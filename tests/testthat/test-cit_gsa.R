# Regression tests for cit_gsa: the asymptotic path plus the permutation path
# that R CMD check's examples did not previously exercise.

make_data <- function(n = 60, r = 8, seed = 1) {
  set.seed(seed)
  X <- data.frame(X = as.factor(rbinom(n, size = 1, prob = 0.5)))
  M <- matrix(rnorm(n * r), nrow = n,
    dimnames = list(NULL, paste0("g", seq_len(r))))
  # make g1 genuinely associated with X (so it survives adaptive thresholds)
  M[, 1] <- M[, 1] + 2 * (as.numeric(X$X) - 1)
  list(M = M, X = X, geneset = list(set1 = paste0("g", 1:4),
    set2 = paste0("g", 5:8)))
}

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

  # permutation: 2 columns only (no test_statistic)
  res_p <- cit_gsa(M = d$M, X = d$X, geneset = d$geneset,
    test = "permutation", n_perm = 20, parallel = FALSE)
  expect_named(res_p$pvals, c("raw_pval", "adj_pval"))
})

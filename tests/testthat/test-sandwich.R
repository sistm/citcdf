# Sandwich (heteroskedasticity-robust) variance of sqrt(n) * beta_hat^X,
# Supplementary Section 1.5.

brute_sigma <- function(D, X, Z = NULL) {
  # explicit construction, one observation at a time, from lm() residuals
  W <- if (is.null(Z)) model.matrix(~., X) else model.matrix(~., cbind(X, Z))
  n <- nrow(W)
  iX <- 1 + seq_len(ncol(model.matrix(~., X)) - 1)
  G <- (solve(crossprod(W) / n) %*% t(W))[iX, , drop = FALSE]   # K x n, col i = gamma_i
  E <- residuals(lm(D ~ W - 1))
  U <- t(vapply(seq_len(n), function(i) as.vector(outer(G[, i], E[i, ])),
    numeric(nrow(G) * ncol(E))))
  crossprod(U) / n
}

test_that(".cit_sandwich_ev() matches a brute-force construction of Sigma_hat", {
  set.seed(11)
  n <- 150
  z <- rnorm(n)
  X <- data.frame(X1 = z + rnorm(n), X2 = factor(sample(letters[1:3], n, TRUE)))
  Z <- data.frame(Z = z)
  Y <- cbind(z + rnorm(n), z + rnorm(n))
  D <- cbind(outer(Y[, 1], c(-1, 0, 1), "<=") * 1, outer(Y[, 2], c(-.5, .5), "<=") * 1)
  for (zz in list(NULL, Z)) {
    design <- .cit_design(X, zz, n)
    ref <- eigen(brute_sigma(D, X, zz), symmetric = TRUE, only.values = TRUE)$values
    ev <- .cit_sandwich_ev(D, design)
    expect_equal(ev, ref[seq_along(ev)], tolerance = 1e-8)
    expect_true(all(abs(ref[-seq_along(ev)]) < 1e-8 * max(ref)))
  }
})

test_that("the n x n Gram path gives the same non-zero spectrum", {
  set.seed(12)
  n <- 30
  X <- data.frame(X = rnorm(n)); Z <- data.frame(Z = rnorm(n))
  D <- outer(rnorm(n), sort(rnorm(40)), "<=") * 1          # 40 columns > n rows
  design <- .cit_design(X, Z, n)
  ev <- .cit_sandwich_ev(D, design)
  E <- qr.resid(design$qr, D); U <- as.vector(design$H) * E
  ref <- eigen(crossprod(U) / n, symmetric = TRUE, only.values = TRUE)$values
  expect_equal(ev, ref[seq_along(ev)], tolerance = 1e-8)
})

test_that("asymptotic tests hold their level when Z affects Y (heteroskedastic errors)", {
  skip_on_cran()
  set.seed(42)
  R <- 400; n <- 200
  p <- replicate(R, {
    z <- rnorm(n); sh <- rnorm(n)
    M <- sapply(1:3, function(g) 1.5 * z + 0.7 * sh + rnorm(n))
    colnames(M) <- paste0("g", 1:3)
    X <- data.frame(X = z + rnorm(n)); Z <- data.frame(Z = z)
    c(cit_asymp(M[, 1], X, Z, space_y = TRUE, number_y = 10)$raw_pval,
      cit_gsa(M, X, Z, geneset = colnames(M), parallel = FALSE)$pvals$raw_pval)
  })
  lvl <- rowMeans(p < 0.05)
  # the former factorised variance gave ~0.005 here (see NEWS)
  expect_true(all(lvl > 0.025 & lvl < 0.08))
})

test_that("residuals = 'restricted' uses the null-model residuals", {
  set.seed(13)
  n <- 120
  z <- rnorm(n)
  X <- data.frame(X = factor(sample(c("a", "b", "c"), n, TRUE)))
  Y <- z + (X$X == "b") + rnorm(n)
  D <- outer(Y, c(-1, 0, 1, 2), "<=") * 1
  for (Z in list(NULL, data.frame(Z = z))) {
    design <- .cit_design(X, Z, n)
    W0 <- if (is.null(Z)) matrix(1, n) else model.matrix(~., Z)
    E0 <- residuals(lm(D ~ W0 - 1))
    G <- t(design$H)
    U <- do.call(cbind, lapply(seq_len(ncol(G)), function(k) G[, k] * E0))
    ref <- eigen(crossprod(U) / n, symmetric = TRUE, only.values = TRUE)$values
    ev <- .cit_sandwich_ev(D, design, residuals = "restricted")
    expect_equal(ev, ref[seq_along(ev)], tolerance = 1e-8)
    # under this alternative, restricted residuals inflate the variance
    expect_gt(sum(ev), sum(.cit_sandwich_ev(D, design, residuals = "full")))
  }
})

test_that("the residuals argument is passed through and leaves the statistic unchanged", {
  set.seed(14)
  n <- 100
  X <- data.frame(X = rnorm(n)); Z <- data.frame(Z = rnorm(n))
  M <- matrix(rnorm(n * 3), n, dimnames = list(NULL, paste0("g", 1:3)))
  a_f <- cit_asymp(M[, 1], X, Z, residuals = "full")
  a_r <- cit_asymp(M[, 1], X, Z, residuals = "restricted")
  expect_equal(a_f$test_statistic, a_r$test_statistic)
  expect_false(isTRUE(all.equal(a_f$raw_pval, a_r$raw_pval)))
  m_r <- cit_multi(M, X, Z, parallel = FALSE, residuals = "restricted")$pvals
  expect_equal(m_r$raw_pval[1],
    cit_asymp(M[, 1], X, Z, space_y = TRUE, residuals = "restricted")$raw_pval)
  g_f <- cit_gsa(M, X, Z, geneset = colnames(M), parallel = FALSE)$pvals
  g_r <- cit_gsa(M, X, Z, geneset = colnames(M), parallel = FALSE,
    residuals = "restricted")$pvals
  expect_equal(g_f$test_statistic, g_r$test_statistic)
  expect_false(isTRUE(all.equal(g_f$raw_pval, g_r$raw_pval)))
  expect_error(cit_asymp(M[, 1], X, Z, residuals = "foo"))
})

test_that("small_sample_corr = \"HC1\" and \"HC2\" rescale the sandwich estimator", {
  set.seed(15)
  n <- 50
  X <- data.frame(X = factor(sample(c("a", "b", "c"), n, TRUE)))
  Z <- data.frame(Z1 = rnorm(n), Z2 = rnorm(n))
  D <- outer(rnorm(n), c(-1, 0, 1), "<=") * 1
  design <- .cit_design(X, Z, n)
  # d = 5 for the full model (1, X, Z), 3 for the null model (1, Z)
  expect_equal(.cit_sandwich_ev(D, design, "full", small_sample_corr = "HC1"),
    .cit_sandwich_ev(D, design, "full") * n / (n - 5))
  expect_equal(.cit_sandwich_ev(D, design, "restricted", small_sample_corr = "HC1"),
    .cit_sandwich_ev(D, design, "restricted") * n / (n - 3))

  # HC2: brute-force Sigma_hat with residuals divided by sqrt(1 - leverage)
  W <- design$modelmat
  G <- t(design$H)
  for (res in c("full", "restricted")) {
    W_res <- switch(res, full = W, restricted = W[, -design$indexes_X])
    P <- W_res %*% solve(crossprod(W_res), t(W_res))
    E <- (D - P %*% D) / sqrt(1 - diag(P))
    U <- do.call(cbind, lapply(seq_len(ncol(G)), FUN = function(k) {
      G[, k] * E
    }))
    ev_ref <- eigen(crossprod(U) / n, symmetric = TRUE, only.values = TRUE)$values
    ev <- .cit_sandwich_ev(D, design, res, small_sample_corr = "HC2")
    expect_equal(ev, ev_ref[seq_along(ev)])
  }
  expect_equal(design$lev, unname(diag(W %*% solve(crossprod(W), t(W)))))
})

test_that("small_sample_corr defaults and is passed through", {
  set.seed(16)
  n <- 100
  X <- data.frame(X = rnorm(n)); Z <- data.frame(Z = rnorm(n))
  M <- matrix(rnorm(n * 3), n, dimnames = list(NULL, paste0("g", 1:3)))
  a <- lapply(c(none = "none", HC1 = "HC1", HC2 = "HC2"), FUN = function(corr) {
    cit_asymp(M[, 1], X, Z, space_y = TRUE, small_sample_corr = corr)
  })
  # "HC2" is the default for a single outcome
  expect_identical(cit_asymp(M[, 1], X, Z, space_y = TRUE), a$HC2)
  expect_equal(a$none$test_statistic, a$HC2$test_statistic)
  expect_gt(a$HC1$raw_pval, a$none$raw_pval)
  expect_gt(a$HC2$raw_pval, a$none$raw_pval)
  expect_false(isTRUE(all.equal(a$HC1$raw_pval, a$HC2$raw_pval)))
  expect_equal(cit_multi(M, X, Z, parallel = FALSE)$pvals$raw_pval[1],
    a$HC2$raw_pval)
  expect_equal(cit_multi(M, X, Z, parallel = FALSE,
    small_sample_corr = "HC1")$pvals$raw_pval[1], a$HC1$raw_pval)
  # "none" is the default for gene sets
  g <- lapply(c(none = "none", HC1 = "HC1", HC2 = "HC2"), FUN = function(corr) {
    suppressMessages(cit_gsa(M, X, Z, geneset = colnames(M), parallel = FALSE,
      small_sample_corr = corr)$pvals)
  })
  expect_identical(cit_gsa(M, X, Z, geneset = colnames(M), parallel = FALSE)$pvals,
    g$none)
  expect_gt(g$HC1$raw_pval, g$none$raw_pval)
  expect_gt(g$HC2$raw_pval, g$none$raw_pval)
  expect_error(cit_asymp(M[, 1], X, Z, small_sample_corr = "yes"))
  expect_error(cit_asymp(M[, 1], X, Z, small_sample_corr = TRUE))
})

test_that("small-sample warnings and message are issued when relevant", {
  set.seed(17)
  n <- 40
  X <- data.frame(X = rnorm(n)); Z <- data.frame(Z = rnorm(n))
  M <- matrix(rnorm(n * 3), n, dimnames = list(NULL, paste0("g", 1:3)))
  n_warn <- function(expr) {
    sum(grepl("anti-conservative", testthat::capture_warnings(expr)))
  }
  # single outcome: d/n = 3/40 > 0.05
  expect_equal(n_warn(cit_asymp(M[, 1], X, Z, small_sample_corr = "none")), 1)
  expect_equal(n_warn(cit_multi(M, X, Z, parallel = FALSE,
    small_sample_corr = "none")), 1)
  expect_no_warning(cit_asymp(M[, 1], X, Z))
  expect_no_warning(cit_asymp(M[, 1], X, Z, small_sample_corr = "HC1"))
  expect_no_warning(cit_asymp(M[, 1], X, Z, residuals = "restricted",
    small_sample_corr = "none"))
  expect_warning(cit_asymp(M[1:20, 1], X[1:20, , drop = FALSE],
    Z[1:20, , drop = FALSE]), "permutation test")
  # gene set: d/n = 3/40 < 0.1, then 5/40 > 0.1
  expect_no_warning(cit_gsa(M, X, Z, geneset = colnames(M), parallel = FALSE))
  expect_message(cit_gsa(M, X, Z, geneset = colnames(M), parallel = FALSE,
    small_sample_corr = "HC1"), "reduces power")
  Z3 <- data.frame(Z1 = rnorm(n), Z2 = rnorm(n), Z3 = rnorm(n))
  expect_warning(cit_gsa(M, X, Z3, geneset = colnames(M), parallel = FALSE),
    "anti-conservative")
  expect_no_message(cit_gsa(M, X, Z3, geneset = colnames(M), parallel = FALSE,
    small_sample_corr = "HC2"))
  # restricted residuals: d/n = 4/20 > 0.1
  expect_equal(n_warn(cit_asymp(M[1:20, 1], X[1:20, , drop = FALSE],
    Z3[1:20, ], residuals = "restricted", small_sample_corr = "none")), 1)
  # the permutation test is not suggested with many covariates
  w <- testthat::capture_warnings(cit_asymp(M[1:20, 1], X[1:20, , drop = FALSE],
    Z3[1:20, ]))
  expect_true(grepl("n = 20 < 30", w) && !grepl("permutation", w))
})

test_that("a warning is issued when a level of X has few observations", {
  set.seed(18)
  n <- 100
  Y <- rnorm(n)
  X <- data.frame(X = factor(rep(c("a", "b"), c(5, n - 5))))
  expect_warning(cit_asymp(Y, X), "few observations \\(here 5\\)")
  expect_no_warning(cit_asymp(Y, X, residuals = "restricted"))
  expect_no_warning(cit_asymp(Y, data.frame(X = factor(rep(c("a", "b"), n / 2)))))
})

test_that("a leverage of 1 gets the HC1 weight under HC2, with a warning", {
  set.seed(19)
  n <- 60
  Y <- rnorm(n)
  X <- data.frame(X = factor(rep(c("a", "b", "c"), c(1, 29, 30))))
  design <- .cit_design(X, NULL, n)
  expect_equal(.cit_hc2_weights(design$lev, n, 3),
    c(n / (n - 3), rep(c(29 / 28, 30 / 29), c(29, 30))))
  w <- testthat::capture_warnings(p_hc2 <- cit_asymp(Y, X)$raw_pval)
  expect_true(any(grepl("1 observation\\(s\\) with a leverage of 1", w)))
  expect_true(is.finite(p_hc2))
  w <- testthat::capture_warnings(cit_asymp(Y, X, small_sample_corr = "HC1"))
  expect_false(any(grepl("leverage", w)))
})

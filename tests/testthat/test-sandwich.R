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

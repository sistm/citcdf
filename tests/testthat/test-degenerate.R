test_that("the permutation test accepts a multi-column X when Z is NULL", {
  set.seed(1)
  n <- 60
  X2 <- data.frame(X1 = as.factor(rbinom(n, 1, 0.5)),
    X2 = rnorm(n))
  M <- data.frame(replicate(4, rnorm(n)))

  # cit_perm() has always accepted this; cit_multi() used to refuse it.
  res <- cit_multi(M, X2, test = "permutation", n_perm = 20,
    adaptive = FALSE, parallel = FALSE)
  expect_s3_class(res, "cit_multi")
  expect_equal(nrow(res$pvals), ncol(M))
  expect_true(all(res$pvals$raw_pval >= 0 & res$pvals$raw_pval <= 1))

  expect_no_error(cit_multi(M, X2, test = "permutation", n_perm = 20,
    adaptive = TRUE, parallel = FALSE))

  # but the conditional scheme is still single-column only
  Z <- data.frame(Z = rnorm(n))
  expect_error(
    cit_multi(M, X2, Z, test = "permutation", n_perm = 10,
      adaptive = FALSE, parallel = FALSE),
    "single column in 'X'")
})

test_that("a constant outcome is rejected up front, not inside eigen()", {
  set.seed(1)
  n <- 60
  X <- data.frame(X = as.factor(rbinom(n, 1, 0.5)))

  for (Y in list(rep(3, n), rep(0, n))) {
    for (sy in c(FALSE, TRUE)) {
      expect_error(cit_asymp(Y, X, space_y = sy), "one single value")
    }
    expect_error(
      cit_perm(Y, X, NULL, X_star = X_perm(X, NULL, n_perm = 10), n_perm = 10),
      "one single value")
  }
})

test_that("cit_multi() and cit_gsa() name the zero-variance outcomes", {
  set.seed(1)
  n <- 60
  X <- data.frame(X = as.factor(rbinom(n, 1, 0.5)))
  M <- data.frame(replicate(4, rnorm(n)))
  colnames(M) <- paste0("g", 1:4)
  M$g2 <- 4
  M$g4 <- 0
  geneset <- list(s1 = c("g1", "g2"), s2 = c("g3", "g4"))

  expect_error(cit_multi(M, X, parallel = FALSE), "g2, g4")
  expect_error(cit_gsa(M = M, X = X, geneset = geneset, parallel = FALSE),
    "g2, g4")
  # matrix input takes the same path
  expect_error(
    cit_gsa(M = as.matrix(M), X = X, geneset = geneset, parallel = FALSE),
    "g2, g4")

  # and a well-behaved M is untouched
  M2 <- data.frame(replicate(4, rnorm(n)))
  expect_no_error(cit_multi(M2, X, parallel = FALSE))
})

test_that("plot.cit_multi() labels the BH line correctly and honours nominal_level", {
  set.seed(1)
  n <- 100
  X <- data.frame(X = as.factor(rbinom(n, 1, 0.5)))
  M <- data.frame(replicate(30, rnorm(n)))
  res <- cit_multi(M, X, test = "asymptotic", parallel = FALSE)

  for (lev in c(0.05, 0.01)) {
    p <- plot(res, nominal_level = lev)
    b <- ggplot2::ggplot_build(p)
    sc <- Filter(function(s) "colour" %in% s$aesthetics, b$plot$scales$scales)[[1]]
    expect_true(any(grepl("Benjamini-Hochberg", sc$labels)))
    expect_false(any(grepl("Bonferroni", sc$labels)))
    # the flat reference line must sit at nominal_level, on the log10 scale
    flat <- Filter(function(d) length(unique(d$y)) == 1L && all(is.finite(d$y)),
      b$data)
    expect_true(any(vapply(flat,
      function(d) isTRUE(all.equal(10^unique(d$y), lev)), logical(1))))
  }
})

test_that("plot_compare_ccdf() legends use the caller's column names", {
  legend_keys <- function(p) {
    out <- character(0)
    walk <- function(q) {
      if (inherits(q, "ggplot")) {
        sc <- Filter(function(s) "colour" %in% s$aesthetics, q$scales$scales)
        if (length(sc)) out <<- c(out, sc[[1]]$labels)
        if (!is.null(q$labels$title)) out <<- c(out, q$labels$title)
      }
      for (r in q$patches$plots) walk(r)
    }
    walk(p)
    out
  }

  set.seed(123)
  n <- 40
  Y  <- data.frame(expression = rnorm(n))
  Xf <- data.frame(treatment = as.factor(rbinom(n, 1, 0.5)))
  Xc <- data.frame(dose = rnorm(n))
  Zf <- data.frame(batch = as.factor(rbinom(n, 1, 0.5)))
  Zc <- data.frame(age = rnorm(n))

  cases <- list(
    list(Y, Xf, NULL, FALSE), list(Y, Xc, NULL, TRUE),
    list(Y, Xf, Zf, FALSE),   list(Y, Xf, Zf, TRUE),
    list(Y, Xc, Zf, FALSE),   list(Y, Xf, Zc, FALSE),
    list(Y, Xc, Zc, FALSE),   list(Y, Xc, Zc, TRUE))
  for (k in seq_along(cases)) {
    a <- cases[[k]]
    keys <- legend_keys(plot_compare_ccdf(a[[1]], a[[2]], a[[3]],
      discretize = a[[4]]))
    info <- paste("case", k)
    # the real names appear, and no bare "X="/"Z=" placeholder survives
    expect_true(any(grepl("treatment|dose", keys)), info = info)
    expect_false(any(grepl("\\bX=|\\bZ=|on X$|on Z$|Given X ", keys)), info = info)
    if (!is.null(a[[3]])) {
      expect_true(any(grepl("batch|age", keys)), info = info)
    }
  }

  # fallback when the column has no usable name
  keys <- legend_keys(plot_compare_ccdf(Y, stats::setNames(Xf, ""),
    discretize = FALSE))
  expect_true(any(grepl("X=", keys)))

  # only the first column name is used, never the whole vector
  p <- plot_compare_ccdf(data.frame(g1 = rnorm(n), g2 = rnorm(n)), Xf,
    discretize = FALSE)
  expect_length(p$labels$x, 1L)
  expect_identical(p$labels$x, "g1")
})

test_that("cit_perm() draws X_star itself when it is not supplied", {
  set.seed(123)
  n <- 100
  X <- data.frame(X = as.factor(rbinom(n, 1, 0.5)))
  Y <- (X$X == 1) * rnorm(n) + (X$X == 0) * rnorm(n, mean = 0.5)
  Zc <- data.frame(Z = rnorm(n))
  Zf <- data.frame(Z = as.factor(rbinom(n, 1, 0.5)))

  # the implicit draw is the same draw, in the same place in the RNG stream
  for (zz in list(NULL, Zc, Zf)) {
    set.seed(7)
    implicit <- cit_perm(Y, X, zz, n_perm = 30)
    set.seed(7)
    explicit <- cit_perm(Y, X, zz, X_star = X_perm(X, zz, n_perm = 30),
      n_perm = 30)
    expect_identical(implicit, explicit)
  }

  # a supplied X_star is used untouched, and a longer pool is truncated
  set.seed(21)
  pool <- X_perm(X, NULL, n_perm = 200)
  expect_identical(cit_perm(Y, X, X_star = pool, n_perm = 50),
    cit_perm(Y, X, X_star = pool[seq_len(50)], n_perm = 50))

  # too short a pool is now a message rather than "subscript out of bounds"
  expect_error(cit_perm(Y, X, X_star = X_perm(X, NULL, n_perm = 5), n_perm = 50),
    "holds 5 permuted design")

  # positional calls written against the old signature still work
  set.seed(7)
  xs <- X_perm(X, NULL, n_perm = 50)
  expect_identical(cit_perm(Y, X, NULL, xs, 50),
    cit_perm(Y, X, NULL, X_star = xs, n_perm = 50))

  # multi-column X without Z goes through X_perm() unchanged
  X2 <- data.frame(X1 = as.factor(rbinom(n, 1, 0.5)), X2 = rnorm(n))
  expect_no_error(cit_perm(Y, X2, n_perm = 20))
})

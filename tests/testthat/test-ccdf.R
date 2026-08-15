test_that("ccdf returns a well-formed object for every method and Z combination", {
  set.seed(1)
  n <- 60
  X <- data.frame(X = as.factor(rbinom(n, 1, 0.5)))
  Z <- data.frame(Z = rnorm(n))
  Y <- rnorm(n)
  for (method in c("OLS", "logistic")) {
    res <- ccdf(Y = Y, X = X, method = method)
    expect_s3_class(res, "ccdf")
    expect_named(res, c("cdf", "ccdf", "y", "x"))
    expect_length(res$ccdf, length(res$y))
    expect_false(anyNA(res$ccdf))

    res_z <- ccdf(Y = Y, X = X, Z = Z, method = method)
    expect_named(res_z, c("cdf", "ccdf_nox", "ccdf_x", "y", "x", "z"))
    expect_false(anyNA(c(res_z$ccdf_x, res_z$ccdf_nox)))
  }
  # the logistic branch returns probabilities
  p <- ccdf(Y = Y, X = X, method = "logistic")$ccdf
  expect_true(all(p >= 0 & p <= 1))
})

test_that("space_y = TRUE uses at most number_y thresholds", {
  # $y holds the observations falling in each threshold bin, not the grid; the
  # grid shows through $cdf, which takes one distinct value per bin.
  set.seed(2)
  n <- 80
  X <- data.frame(X = as.factor(rbinom(n, 1, 0.5)))
  Y <- rnorm(n)
  full <- ccdf(Y = Y, X = X, space_y = FALSE)
  grid <- ccdf(Y = Y, X = X, space_y = TRUE, number_y = 10)
  expect_equal(length(unique(full$cdf)), length(unique(Y)) - 1L)
  expect_lte(length(unique(grid$cdf)), 9L)
  expect_lte(max(grid$y), max(Y))
})

test_that("ccdf handles degenerate inputs and ignores the row order", {
  set.seed(3)
  n <- 40
  X <- data.frame(X = as.factor(rbinom(n, 1, 0.5)))
  Z <- data.frame(Z = rnorm(n))
  expect_no_error(ccdf(Y = rep(2, n), X = X, space_y = TRUE, number_y = 5))
  single <- ccdf(Y = seq_len(n) + rnorm(n, sd = 1e-6), X = X, space_y = FALSE)
  expect_length(single$ccdf, n - 1L)
  expect_false(anyNA(single$ccdf))

  # permuting the observations permutes the output, nothing more
  Y <- rnorm(n)
  o <- sample.int(n)
  a <- ccdf(Y = Y, X = X, Z = Z)
  b <- ccdf(Y = Y[o], X = X[o, , drop = FALSE], Z = Z[o, , drop = FALSE])
  expect_equal(a$ccdf_x[order(a$y)], b$ccdf_x[order(b$y)])
})

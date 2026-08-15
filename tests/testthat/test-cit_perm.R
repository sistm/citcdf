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
  # score is a count bounded by n_perm (not the statistic)
  expect_true(res$score >= 0 && res$score <= 30)
})

test_that("cit_asymp and cit_perm report the SAME observed statistic", {
  set.seed(5)
  n <- 100
  X <- data.frame(X = as.factor(rbinom(n, 1, .5)))
  Y <- rnorm(n) + 0.6 * (as.numeric(X$X) - 1)
  Xs <- X_perm(X, NULL, n_perm = 30)

  # space_y = FALSE so both use sort(unique(Y))
  a <- cit_asymp(Y, X, space_y = FALSE)
  p <- cit_perm(Y, X, NULL, X_star = Xs, n_perm = 30, space_y = FALSE)

  expect_identical(a$test_statistic, p$test_statistic)
})

test_that("X_perm preserves X's type: class and levels for a factor, values for a numeric", {
  set.seed(20)
  n <- 60
  Zs <- list(cont = data.frame(Z = rnorm(n)),
    disc = data.frame(Z = factor(sample(c("a", "b"), n, TRUE))),
    none = NULL)
  Xs <- list(labelled = factor(ifelse(rbinom(n, 1, 0.5) == 1, "treated", "control")),
    numlevel = factor(rbinom(n, 1, 0.5)),
    numeric  = rbinom(n, 1, 0.5))
  for (xn in names(Xs)) for (zn in names(Zs)) {
    lab <- paste(xn, zn)
    X <- Xs[[xn]]
    xs <- X_perm(data.frame(X = X), Zs[[zn]], n_perm = 3)[[1]]
    xs <- if (is.data.frame(xs)) xs[, 1] else xs
    expect_false(anyNA(xs), info = lab)
    if (is.factor(X)) {
      expect_s3_class(xs, "factor")
      expect_identical(levels(xs), levels(X), info = lab)
    } else {
      expect_true(is.numeric(xs), info = lab)
    }
    # a permutation preserves the marginal exactly
    expect_identical(sort(as.character(xs)), sort(as.character(X)), info = lab)
  }
})

test_that("X_perm dispatches on the kind of Z and preserves discrete strata", {
  set.seed(30)
  n <- 90
  X <- data.frame(X = as.factor(rbinom(n, 1, 0.5)))
  Zs <- list(numeric = rnorm(n),
    double  = as.double(rnorm(n)),
    integer = sample.int(5, n, TRUE),
    logical = rbinom(n, 1, 0.5) == 1,
    factor  = factor(sample(c("a", "b", "c"), n, TRUE)),
    classed = structure(rnorm(n), class = c("labelled", "numeric")))
  for (nm in names(Zs)) {
    z <- Zs[[nm]]
    xs <- X_perm(X, data.frame(Z = z), n_perm = 2)
    expect_length(xs, 2)
    expect_identical(sort(as.character(xs[[1]][, 1])), sort(as.character(X$X)),
      info = nm)
    if (!is.numeric(z)) {
      # permuting within strata leaves every stratum's composition untouched
      expect_identical(table(xs[[1]][, 1], z, dnn = NULL),
        table(X$X, z, dnn = NULL), info = nm)
    }
  }
  expect_error(X_perm(X, data.frame(Z = I(as.list(seq_len(n)))), n_perm = 1),
    "unsupported class")
})

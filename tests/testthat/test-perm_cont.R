test_that("perm_cont returns a permutation of X, not a resample", {
  set.seed(1)
  n <- 60
  Z <- rnorm(n)
  X <- rbinom(n, 1, plogis(0.8 * Z))
  # a permutation preserves the whole marginal distribution exactly, every draw
  expect_true(all(replicate(20, identical(sort(perm_cont(X, Z)), sort(X)))))
})

test_that("perm_cont draws a bijection of the donor indices", {
  set.seed(2)
  n <- 50
  Z <- rnorm(n)
  X <- rbinom(n, 1, plogis(0.8 * Z))
  # feeding the indices themselves back through the sampler exposes the donors
  expect_true(all(replicate(20, identical(sort(perm_cont(seq_len(n), Z)), seq_len(n)))))
})

test_that("perm_cont still matches donors on Z", {
  set.seed(3)
  n <- 200
  Z <- rnorm(n)
  X <- rbinom(n, 1, plogis(0.8 * Z))
  mm <- model.matrix(~Z)
  fit <- as.vector(mm %*% (solve(crossprod(mm)) %*% t(mm) %*% X))
  matched <- mean(replicate(20, mean(abs(fit - fit[perm_cont(seq_len(n), Z)]))))
  uncond  <- mean(replicate(20, mean(abs(fit - fit[sample.int(n)]))))
  # conditioning on Z should leave donors far better matched than a plain
  # permutation; the measured ratio is ~0.05, so 0.5 is a loose guard
  expect_lt(matched, 0.5 * uncond)
})

test_that("perm_cont never returns a degenerate (single-level) design", {
  # Drawing donors WITH replacement can return an X_star with only one level,
  # which then kills cit_perm() in model.matrix() ("contrasts can be applied
  # only to factors with 2 or more levels"). Measured at n = 20 with a Gaussian
  # kernel and replacement: 0.05% of draws, i.e. 23 of 600 tests failing over
  # n_perm = 299 draws each. Drawing without replacement makes it impossible.
  set.seed(5)
  ok <- vapply(c(20L, 30L), function(n) {
    Z <- rnorm(n)
    X <- rbinom(n, 1, plogis(0.8 * Z))
    if (length(unique(X)) < 2L) return(TRUE)
    all(replicate(20, length(unique(perm_cont(X, Z))) == length(unique(X))))
  }, TRUE)
  expect_true(all(ok))
})

test_that("perm_cont copes with a degenerate fit", {
  set.seed(4)
  n <- 30
  Z <- rnorm(n)
  X <- rep(1, n)          # fitted values all identical -> all distances zero
  xs <- perm_cont(X, Z)
  expect_identical(sort(xs), sort(X))
})

test_that("perm_cont keeps conditioning when the Gaussian weights underflow", {
  # exp(-d^2/(2h^2)) underflows to 0 once d/h > ~38.6. A leverage point in Z
  # puts every donor beyond that, and a uniform fallback would hand exactly the
  # most influential observation a completely unconditioned draw. It must take
  # the nearest donor instead.
  set.seed(6)
  n <- 200
  Z <- c(rnorm(n - 1), 40)                       # one far-out covariate value
  X <- rbinom(n, 1, plogis(0.8 * pmin(Z, 3)))
  mm <- model.matrix(~Z)
  fit <- as.vector(mm %*% (solve(crossprod(mm)) %*% t(mm) %*% X))
  h <- stats::sd(fit) * n^(-1/3)
  expect_true(all(exp(-(fit[n] - fit[-n])^2 / (2 * h^2)) == 0))   # really underflowed
  nearest <- order(abs(fit[n] - fit))[2]                          # [1] is itself
  donors <- replicate(30, perm_cont(seq_len(n), Z)[n])
  # a uniform draw would land on the nearest donor ~1/199 of the time
  expect_true(mean(donors == nearest) > 0.5)
})

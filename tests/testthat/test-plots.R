expect_buildable <- function(p) {
  expect_true(inherits(p, "ggplot") || inherits(p, "patchwork"))
  expect_no_error(invisible(ggplot2::ggplot_build(p)))
}

test_that("plot_compare_ccdf builds for every X/Z combination it dispatches on", {
  skip_if_not_installed("ggplot2")
  set.seed(1)
  n <- 40
  Y <- data.frame(Y = rnorm(n))
  Xf <- data.frame(X = as.factor(rbinom(n, 1, 0.5)))
  Xc <- data.frame(X = rnorm(n))
  X3 <- data.frame(X = as.factor(sample(0:2, n, TRUE)))
  Zf <- data.frame(Z = as.factor(rbinom(n, 1, 0.5)))
  Zc <- data.frame(Z = rnorm(n))
  for (cs in list(list(Xf, NULL), list(Xf, Zf), list(Xf, Zc),
    list(Xc, Zc), list(X3, NULL))) {
    expect_buildable(plot_compare_ccdf(Y, cs[[1]], cs[[2]]))
  }
})

test_that("the plot methods for cit_multi and cit_gsa build", {
  skip_if_not_installed("ggplot2")
  d <- make_data(n = 50, r = 8)
  expect_buildable(plot(cit_multi(M = data.frame(d$M), X = d$X,
    test = "asymptotic", parallel = FALSE)))
  gsa <- cit_gsa(M = data.frame(d$M), X = d$X, geneset = d$geneset,
    test = "asymptotic", parallel = FALSE)
  expect_buildable(plot(gsa, M = data.frame(d$M), X = d$X,
    geneset = d$geneset))
})

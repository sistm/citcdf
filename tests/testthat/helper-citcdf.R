# shared fixtures for the citcdf test suite
make_data <- function(n = 60, r = 8, seed = 1) {
  set.seed(seed)
  X <- data.frame(X = as.factor(rbinom(n, size = 1, prob = 0.5)))
  M <- matrix(rnorm(n * r), nrow = n,
    dimnames = list(NULL, paste0("g", seq_len(r))))
  M[, 1] <- M[, 1] + 2 * (as.numeric(X$X) - 1)
  list(M = M, X = X, geneset = list(set1 = paste0("g", 1:4),
    set2 = paste0("g", 5:8)))
}

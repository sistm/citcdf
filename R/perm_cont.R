#' Permutation procedure when Z is continuous
#'
#' @param X a numeric or factor vector of length \code{n}
#' containing the variable to be tested (the condition to be tested).
#'
#' @param Z a numeric vector of length \code{n}
#' containing the covariate. Multiple variables are not allowed.
#'
#' @seealso \code{\link{X_perm}}
#'
#' @export
#'
#' @import stats
#'
#' @return \code{X_star} a vector of permuted \code{X}.
#'
#' @examples
#'
#' if (interactive()) {
#'   X <- rbinom(n = 100, size = 1, prob = 0.5)
#'   Z <- rnorm(100, 0, 1)
#'   Y <- ((X == 1) * rnorm(n = 50, 0, 1)) + ((X == 0) * rnorm(n = 50, 0.5, 1))
#'   res <- perm_cont(X, Z)
#' }
perm_cont <- function(X, Z) {
  modmat <- model.matrix(~Z)
  reg_coefs <- solve(crossprod(modmat)) %*% t(modmat) %*% X
  X_star <- rep(NA, length(X))
  n <- length(Z)

  fit <- as.vector(modmat %*% reg_coefs)
  for (i in seq_len(n)) {
    # Weights are proportional to 1/(difference in fitted values)^2.
    # If pred_z^2 is constant within a row, it cancels when sample() normalises the
    # weights, but would make the row all-zero when pred_z == 0.
    w <- 1 / (fit[i] - fit[-i])^2
    ties <- !is.finite(w)
    # Observations sharing i's fitted value then carry infinite weights.
    # But actually, in the limit, all the mass sits on them, uniformly.
    # Sample within the tied set only then.
    if (any(ties)) {
      w <- as.numeric(ties)
    }
    X_star[i] <- sample(X[-i], size = 1, prob = w)
  }

  return(X_star)
}

#' Asymptotic test for conditional independence
#'
#' Test the conditional independence of Y and X given Z.
#'
#' @param Y a numeric vector of length \code{n} to test for conditional independence
#' with \code{X} adjusted on \code{Z}
#'
#' @param X a data frame of size \code{n x p} of numeric or factor vector(s)
#' containing the variable(s) to be tested for conditional independence
#' against \code{X} adjusted on \code{Z}.
#'
#' @param Z a data frame of size \code{n x q} of numeric or factor vector(s)
#' containing the covariate(s) to condition the independence
#' test upon.
#'
#' @param space_y a logical flag indicating whether the y thresholds are spaced.
#' When \code{space_y} is \code{TRUE}, a regular sequence between the minimum and
#' the maximum of the observations is used. Default is \code{FALSE}.
#'
#' @param number_y an integer value indicating the number of y thresholds (and therefore
#' the number of regressions) to perform the test. Default is \code{length(unique(Y))}.
#'
#' @param design an optional (and technical) list of design quantities, as returned by the
#' internal \code{.cit_design(X, Z, n)}. This is used by \code{cit_multi()},
#' to loop-call over many genes while building the model matrix, and computing its
#' cross-product and its inverse only once.
#' Default is \code{NULL}, in which case they are computed from \code{X} and
#' \code{Z}. Users should not be using this argument
#'
#' @importFrom survey pchisqsum
#'
#' @seealso \code{\link{cit_perm}}, \code{\link{cit_multi}}, \code{\link{ccdf}}
#'
#' @export
#'
#' @return A data frame with the following elements:
#' \itemize{
#'   \item \code{raw_pval} contains the raw p-values for a given gene.
#'   \item \code{test_statistic} contains the test statistic for a given gene.
#' }
#'
#' @examples
#'
#' X <- as.factor(rbinom(n = 100, size = 1, prob = 0.5))
#' Y <- ((X == 1) * rnorm(n = 100, 0, 1)) + ((X == 0) * rnorm(n = 100, 0.5, 1))
#' res_asymp <- cit_asymp(Y, data.frame(X = X))
#'
#'
#' Z <- as.factor(rbinom(n = 100, size = 1, prob = 0.5))
#' X <- as.numeric(Z) - 1  + rnorm(n = 100, sd = 1)
#' r <- 1000
#' Y <- replicate(r, as.numeric(Z) - 1)
#' YY <- (Y == 1) * rnorm(n = 100 * r, 0, 1) + (Y == 0) * rnorm(n = 100 * r, 0.5, 1)
#' pvals_sim <- pbapply::pbsapply(1:1000, function(i) {
#'   res_asymp <- cit_asymp(YY[, i], data.frame(X = X), data.frame(Z = Z))
#'   return(res_asymp$raw_pval)
#' })
#' hist(pvals_sim)
#' quantile(pvals_sim)
#'
cit_asymp <- function(Y, X, Z = NULL, space_y = FALSE, number_y = length(unique(Y)),
                      design = NULL) {
  # Quantities that depend only on (X, Z), not on Y. Callers looping over many
  # genes (cit_multi) build this once and pass it in; a direct call computes it
  # on demand, so the public behaviour is unchanged.
  n_Y_all <- length(Y)
  stopifnot(nrow(X) == n_Y_all)
  stopifnot(is.null(Z) || nrow(Z) == n_Y_all)
  if (is.null(design)) {
    design <- .cit_design(X, Z, n_Y_all)
  }
  H <- design$H
  # computing the test statistic
  # depends on Y: has to be recomputed for each gene
  Y <- as.numeric(Y) # is this really necessary ??
  oY <- order(Y)
  y <- .cit_y_grid(Y, space_y, number_y)
  p <- length(y) # number of thresholds used

  index_jumps <- findInterval(y[-p], Y[oY])
  beta <- c(apply(X = H[, oY, drop = FALSE], MARGIN = 1, FUN = cumsum)[index_jumps, , drop = FALSE]) / n_Y_all
  test_stat <- sum(beta^2) * n_Y_all

  # Computing the variance ----
  prop <- index_jumps / n_Y_all

  B <- prop - prop %x% t(prop)
  Bsym  <- B * upper.tri(B, diag = TRUE) + t(B * upper.tri(B, diag = FALSE))
  Sigma <- (tcrossprod(H) / n_Y_all) %x% Bsym

  decomp <- eigen(Sigma, symmetric = TRUE, only.values = TRUE)


  # computing the pvalue ----
  pval <- try(survey::pchisqsum(test_stat, lower.tail = FALSE, df = rep(1, ncol(Sigma)),
    a = decomp$values, method = "saddlepoint"),
  silent = TRUE)
  if (inherits(pval, "try-error")) {
    pval <- try(survey::pchisqsum(test_stat, lower.tail = FALSE, df = rep(1, ncol(Sigma)),
      a = decomp$values, method = "satterthwaite"),
    silent = TRUE)
    if (inherits(pval, "try-error")) {
      pval <- NA
    }
  }

  return(data.frame("raw_pval" = pval, "test_statistic" = test_stat))

}

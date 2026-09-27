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
#' the number of regressions) to perform the test. Only used if \code{space_y}
#' is \code{TRUE}. Default is \code{10}.
#'
#' @param design an optional (and technical) list of design quantities, as returned by the
#' internal \code{.cit_design(X, Z, n)}. This is used by \code{cit_multi()},
#' to loop-call over many genes while building the model matrix, and computing its
#' cross-product and its inverse only once.
#' Default is \code{NULL}, in which case they are computed from \code{X} and
#' \code{Z}. Users should not be using this argument
#'
#' @param variance a character string, the estimator of the covariance
#' of the OLS coefficients of \code{X} that gives the weights of the
#' asymptotic \eqn{\chi^2} mixture. Either \describe{
#'   \item{\code{"sandwich"}}{(default) use the heteroskedasticity-robust
#'   sandwich estimator, valid with or without \code{Z}.}
#'   \item{\code{"independent"}}{the factorized
#'   estimator only valid only when \code{Z} is NULL \code{Y} and
#'   \code{Z} are independent ; too conservative when \code{Z} actually affects \code{Y}.}
#' }
#'
#' @param residuals a character string indicating which residuals are used in
#' the sandwich estimator (ignored when \code{variance = "independent"}). Either \describe{
#'   \item{\code{"full"}}{(default) residuals of the linear model including
#'   \code{X} (Wald-type). Consistent for the variance of the \code{X}
#'   coefficients under the null hypothesis and under the alternative.}
#'   \item{\code{"restricted"}}{residuals of the null model without
#'   \code{X} (score-type). Consistent under the null hypothesis only: under
#'   the alternative the effect of \code{X} is counted as noise, which can
#'   decrease power. Can be quite conservative in small samples, especially if
#'   a level of a factor from \code{X} contains few observations. Can serve as
#'   a conservative sensitivity analysis.}
#' }
#' Only used by the asymptotic test.
#'
#' @importFrom survey pchisqsum
#'
#' @details The \code{space_y} / \code{number_y} grid controls both the
#' resolution of the statistic and its computational cost. See
#' \code{\link{cit_multi}} for details on this trade-off.
#'
#'
#' @return A data frame with the following elements:
#' \itemize{
#'   \item \code{raw_pval} contains the raw p-values for a given gene.
#'   \item \code{test_statistic} contains the test statistic for a given gene.
#' }
#'
#' @seealso \code{\link{cit_perm}}, \code{\link{cit_multi}}, \code{\link{ccdf}}
#'
#' @references Gauthier M, Agniel D, Thiébaut R & Hejblum BP (2021).
#' Distribution-free complex hypothesis testing for single-cell RNA-seq
#' differential expression analysis, \emph{bioRxiv} 445165.
#' \doi{10.1101/2021.05.21.445165}.
#'
#' @export
#'
#' @examples
#'
#' set.seed(123)
#' X <- as.factor(rbinom(n = 100, size = 1, prob = 0.5))
#' Y <- ((X == 1) * rnorm(n = 100, 0, 1)) + ((X == 0) * rnorm(n = 100, 0.5, 1))
#' res_asymp <- cit_asymp(Y, data.frame(X = X))
#'
#'
#' Z <- as.factor(rbinom(n = 100, size = 1, prob = 0.5))
#' X <- as.numeric(Z) - 1  + rnorm(n = 100, sd = 1)
#' r <- 500
#' Y <- replicate(r, as.numeric(Z) - 1)
#' YY <- (Y == 1) * rnorm(n = 100 * r, 0, 1) + (Y == 0) * rnorm(n = 100 * r, 0.5, 1)
#' pvals_sim <- sapply(seq_len(r), function(i) {
#'   cit_asymp(YY[, i], data.frame(X = X), data.frame(Z = Z))$raw_pval
#' })
#' hist(pvals_sim) # well calibrated p-values are uniform under the null
#' quantile(pvals_sim)
#'
cit_asymp <- function(Y, X, Z = NULL, space_y = FALSE, number_y = 10,
                      design = NULL, variance = c("sandwich", "independent"),
                      residuals = c("full", "restricted")) {
  variance <- match.arg(variance)
  residuals <- match.arg(residuals)
  # when 'design' is supplied from `cit_multi`, it has already warned once
  if (variance == "independent" && is.null(design)) {
    .cit_warn_independent()
  }
  # Quantities that depend only on (X, Z), not on Y. Callers looping over many
  # genes (cit_multi) build this once and pass it in; a direct call computes it
  # on the fly.
  n_Y_all <- length(Y)
  stopifnot(nrow(X) == n_Y_all)
  stopifnot(is.null(Z) || nrow(Z) == n_Y_all)
  .cit_check_Y(Y)
  if (is.null(design)) {
    design <- .cit_design(X, Z, n_Y_all)
  }
  H <- design$H


  # threshold indicators D (n x (p-1)): depend on Y, so it needs to be
  # recomputed for each gene.
  Y <- as.numeric(Y) # is this really necessary ?? Or should it be part of .cit_check_Y(Y) ?
  y <- .cit_y_grid(Y, space_y, number_y)
  p <- length(y) # number of thresholds used; the last one (D == 1) is dropped
  D <- outer(Y, y[-p], "<=") * 1

  # computing the test statistic ----
  # beta_hat^X = (W'W)^{-1} W' D restricted to the X rows = H D / n
  beta <- (H %*% D) / n_Y_all
  test_stat <- sum(beta^2) * n_Y_all

  # Computing the eigen values from the empirical variance ----
  ev <- switch(variance,
    sandwich    = .cit_sandwich_ev(D, design, residuals),
    independent = .cit_independent_ev(D, design))

  # computing the pvalue ----
  pval <- try(survey::pchisqsum(test_stat, lower.tail = FALSE, df = rep(1, length(ev)),
    a = ev, method = "saddlepoint"),
  silent = TRUE)
  if (inherits(pval, "try-error")) {
    pval <- try(survey::pchisqsum(test_stat, lower.tail = FALSE, df = rep(1, length(ev)),
      a = ev, method = "satterthwaite"),
    silent = TRUE)
    if (inherits(pval, "try-error")) {
      pval <- NA
    }
  }

  return(data.frame("raw_pval" = pval, "test_statistic" = test_stat))

}

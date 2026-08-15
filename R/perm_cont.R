#' Permutation procedure when Z is continuous
#'
#' @param X a numeric or factor vector of length \code{n}
#' containing the variable to be tested (the condition to be tested).
#'
#' @param Z a numeric vector of length \code{n}
#' containing the covariate. Multiple variables are not allowed.
#'
#' @details
#' \code{X_star} is a conditional permutation of \code{X}: draws are
#' without replacement, so \code{X_star} is a bijective rearrangement of
#' \code{X} and its marginal distribution is preserved exactly. The permut for
#' observation \code{i} is drawn with a weight that decreases with the distance
#' between fitted values of \code{X} given \code{Z}, so they are matched on
#' \code{Z}.
#'
#' For binary \code{X} the fitted value is a linear-probability approximation of the
#' propensity score (matching on the propensity score balances \code{Z},
#' [Rosenbaum and Rubin 1983] and allows to match on that single
#' number rather than on \code{Z}). Of note, this approximation can fall
#' outside \code{[0, 1]} when the association between \code{X} and
#' \code{Z} is strong.
#'
#' This is the conditional permutation test (CPT) of Berrett et al. (2020), which
#' draws a permutation with probability proportional to the likelihood of the permuted
#' assignment under the conditional law of \code{X} given \code{Z}. The
#' sequential draw used here approximates that distribution rather than sampling
#' from it exactly. Permuting rather than resampling \code{X}, as the
#' conditional randomisation test (CRT) of Candes et al. (2018) does, preserves
#' the empirical marginal distribution of \code{X} exactly regardless of the
#' error in the fitted conditional law model (and never produces degenerate
#' single-level designs).
#'
#' Permutation weights are Gaussian in the distance between fitted values, with
#' bandwidth \code{sd(fit) * n^(-1/3)} (polynomial weights leave too much mass on
#' distant candidates, so the neighbourhood was not local and the
#' \code{X}-\code{Z} relationship was not preserved). Of note, the bandwith
#' exponent is larger than Silverman's \code{1/5} because his rate is
#' optimized for density estimation but leaves a first-order bias that
#' invalidates inference built on it (Hall, 1992; Armstrong and Kolesar, 2020);
#' \code{o(n^(-1/4))} is the condition Kim et al. (2022) prove for the
#' analogous local permutation test.
#'
#' @references
#' Berrett TB, Wang Y, Barber RF, Samworth RJ (2020).
#' The conditional permutation test for independence while controlling for
#' confounders. \emph{Journal of the Royal Statistical Society Series B},
#' \bold{82}(1), 175-197.
#'
#' Candes E, Fan Y, Janson L, Lv J (2018). Panning for gold: 'model-X' knockoffs
#' for high dimensional controlled variable selection.
#' \emph{Journal of the Royal Statistical Society Series B}, \bold{80}(3), 551-577.
#'
#' Hemerik J, Goeman JJ (2018). Exact testing with random permutations.
#' \emph{TEST}, \bold{27}(4), 811-825.
#'
#' Rosenbaum PR, Rubin DB (1983). The central role of the propensity score in
#' observational studies for causal effects. \emph{Biometrika}, \bold{70}(1),
#' 41-55.
#'
#' Kim I, Neykov M, Balakrishnan S, Wasserman L (2022).
#' Local permutation tests for conditional independence.
#' \emph{The Annals of Statistics}, \bold{50}(6), 3388-3414.
#'
#' Hall P (1992). Effect of bias estimation on coverage accuracy of bootstrap
#' confidence intervals for a probability density.
#' \emph{The Annals of Statistics}, \bold{20}(2), 675-694.
#'
#' Armstrong TB, Kolesar M (2020). Simple and honest confidence intervals in
#' nonparametric regression. \emph{Quantitative Economics}, \bold{11}(1), 1-39.
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
#' set.seed(123)
#' X <- rbinom(n = 100, size = 1, prob = 0.5)
#' Z <- rnorm(100, 0, 1)
#' X_star <- perm_cont(X, Z)
#' table(X, X_star)
perm_cont <- function(X, Z) {
  modmat <- model.matrix(~Z)
  n <- length(Z)
  fit <- as.vector(X - .lm.fit(modmat, X)$residuals)

  # Data driven bandwidth of the Gaussian kernel.
  h <- stats::sd(fit) * n^(-1 / 3)

  # Draws are WITHOUT replacement, so X_star is a genuine permutation of
  # X rather than a resample (drawing each i independently let 2 observations
  # take the same replacements and then sum(X_star) drifted away from sum(X),
  # leading to non exchangeable designs between the permuted and the observed ones.
  # Of note: the visiting order is randomized (with a fixed order the last
  # observations would systematically get the worst of the remaining candidates...)
  permut <- integer(n)
  available <- rep(TRUE, n)
  ord <- sample.int(n)

  for (pos in seq_along(ord)) {
    i <- ord[pos]
    candidates <- which(available)
    candidates <- candidates[candidates != i]
    if (length(candidates) == 0L) {
      # only i itself is left: permut with an already-assigned position
      j <- ord[sample.int(pos - 1L, 1L)]
      permut[i] <- permut[j]
      permut[j] <- i
      available[i] <- FALSE
      next
    }
    # Weights are Gaussian in the difference in fitted values. The tail is what
    # matters as previous 1/d^2 (even with an added ridge penalty 1/(d^2 + h^2))
    # left too much mass on distant observations and X-Z relationship was not
    # preserved by the permutation.
    d <- abs(fit[i] - fit[candidates])
    w <- exp(-d^2 / (2 * h^2))
    w[!is.finite(w)] <- 0

    j <- if (length(candidates) == 1L) {
      candidates
    } else if (any(w > 0)) {
      candidates[sample.int(length(candidates), 1L, prob = w)]
    } else if (h > 0) {
      # Every candidate is more than ~38 bandwidths away and the weights have all
      # underflowed to 0: neighbourhood is EMPTY (not the same as all candidates
      # being equally good).
      # => take the nearest one.
      # Of note: a uniform draw here would un-condition
      # exactly the highest-leverage observations...
      candidates[which.min(d)]
    } else {
      # h == 0: fitted values are constant (e.g. X is constant), so there
      # really is nothing to condition on and a uniform draw is then correct.
      candidates[sample.int(length(candidates), 1L)]
    }
    permut[i] <- j
    available[j] <- FALSE
  }

  return(X[permut])
}

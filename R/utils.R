#' Shared internal utility functions
#'
#' Definition of the y-threshold grid used by ccdf(), cit_asymp(), cit_perm()
#' and cit_gsa(), so that the asymptotic and permutation tests are evaluate the
#' process on the same thresholds. Callers always drop the last grid point
#' (`y[-p]`), so `to = max(Y)` is the first threshold NOT used: the degenerate
#' threshold at which the indicator is identically 1 is already excluded by the
#' loop.
#
#' @keywords internal
#' @noRd
.cit_y_grid <- function(Y, space_y, number_y) {
  Y <- as.numeric(Y)
  if (!space_y) {
    return(sort(unique(Y)))
  }

  nz <- Y[Y != 0]
  from <- ifelse(length(nz) == 0L, min(Y), min(nz))
  return(seq(from = from, to = max(Y), length.out = number_y))
}

# Degenerate-outcome guard. A Y with only one single value has no conditional
# CDF to estimate, and the failure could be silent or obscure (depending on the
# grid):
#  - space_y = FALSE: .cit_y_grid() returns a length-1 grid, y[-p] is empty,
#    Sigma comes out 0 x 0 and eigen() aborts with "0 x 0 matrix";
#  - space_y = TRUE:  the grid collapses onto a single point, the statistic is
#    ~1e-30 and pchisqsum() returns a silent NaN p-value;
#  - space_y = TRUE and Y identically 0: min(Y[-which(Y == 0)]) is Inf and
#    seq() aborts with "'from' must be a finite number".
# All 3 cases are caught here instead.
.cit_is_constant <- function(u) {
  u <- as.numeric(u)
  # ~3x cheaper than length(unique(u)) < 2 on a 500 x 20000 M, and treats a
  # column whose only observed value is repeated (1, NA, 1) as constant.
  !any(u != u[1L], na.rm = TRUE)
}

.cit_check_Y <- function(Y, arg = "Y") {
  if (.cit_is_constant(Y)) {
    stop("'", arg, "' has one single value. A conditional CDF cannot be ",
      "estimated from a constant outcome, and the test statistic is degenerate. ",
      "Remove zero-variance outcomes before testing.", call. = FALSE)
  }
  invisible(TRUE)
}

# Column-wise version of .cit_check_Y() for the cit_multi() / cit_gsa() inputs.
# A whole matrix is thus only  screened once, plus the offending columns are
# named. M can be a data.frame or a matrix.
.cit_check_M <- function(M) {
  get_col <- if (is.data.frame(M)) function(j) M[[j]] else function(j) M[, j]
  r <- ncol(M)
  const <- vapply(seq_len(r),
    function(j) .cit_is_constant(get_col(j)),
    FUN.VALUE = logical(1))
  if (any(const)) {
    nms <- colnames(M)[const]
    if (is.null(nms)) {
      nms <- paste0("column ", which(const))
    }
    shown <- nms[seq_len(min(5L, length(nms)))]
    stop(sum(const), " of the ", r, " outcomes in 'M' has one single ",
      "value (zero variance): ", paste(shown, collapse = ", "),
      if (length(nms) > 5L) paste0(", ... and ", length(nms) - 5L, " more") else "",
      ".\n  A conditional CDF cannot be estimated from a constant outcome. ",
      "Remove them first, e.g.:\n",
      "  M <- M[, apply(M, 2, function(u) length(unique(u)) > 1L), drop = FALSE]",
      call. = FALSE)
  }
  invisible(TRUE)
}

# Design quantities depending only on (X, Z): constant across genes and,
# when Z is absent, across permutations of X.
.cit_design <- function(X, Z = NULL, n) {
  colnames(X) <- paste0("X", seq_len(ncol(X)))
  if (is.null(Z)) {
    modelmat <- model.matrix(~., data = X)
  } else {
    colnames(Z) <- paste0("Z", seq_len(ncol(Z)))
    modelmat <- model.matrix(~., data = cbind(X, Z))
  }
  indexes_X <- which(substring(colnames(modelmat), 1, 1) == "X")
  H <- n * (solve(crossprod(modelmat)) %*% t(modelmat))[indexes_X, , drop = FALSE]
  # crossprod(modelmat) is invariant under row permutation when Z is absent,
  # so its inverse (restricted to the X rows) is reusable for every permuted design.
  XtXinv_X <- if (is.null(Z)) solve(crossprod(modelmat))[indexes_X, , drop = FALSE] else NULL

  # compute the QR decomposition of the full design (1, X, Z) and of the null
  # design (1, Z) (reused for the OLS residuals of every gene/threshold in the
  # sandwich variance .cit_sandwich_ev(), for residuals = "full" and
  # "restricted" respectively).
  # lev and lev0 are the leverages of the two designs.
  qr_full <- qr(modelmat)
  qr_null <- qr(modelmat[, -indexes_X, drop = FALSE])

  return(list(modelmat = modelmat, indexes_X = indexes_X, H = H,
    XtXinv_X = XtXinv_X, qr = qr_full, qr0 = qr_null,
    lev = rowSums(qr.Q(qr_full)^2), lev0 = rowSums(qr.Q(qr_null)^2))
  )
}

# Eigenvalues of the heteroskedasticity-robust sandwich estimate of the
# asymptotic covariance Sigma, computed from the OLS residuals of the threshold
# indicators D on the full design (1, X, Z) if residuals = "full", or on the
# null design (1, Z) if residuals = "restricted".
# Only the eigenvalues are needed and the test statistic is a squared norm, so
# the ordering of the columns of U is irrelevant.
# When U has more columns than rows, the non-zero eigenvalues are obtained from
# the n x n Gram matrix U U^T / n instead (same non-zero spectrum, cheaper).
# small_sample_corr = "HC1" multiplies Sigma_hat by n / (n - d), with d the
# number of coefficients of the model used when computing the residuals.
# small_sample_corr = "HC2" divides each row of the residuals by
# sqrt(1 - h_i) instead (h_i being the observation levergae).
#
# D: n x m matrix of threshold indicators (m = thresholds x genes).
.cit_sandwich_ev <- function(D, design, residuals = c("full", "restricted"),
                             small_sample_corr = c("none", "HC1", "HC2")) {
  residuals <- match.arg(residuals)
  small_sample_corr <- match.arg(small_sample_corr)
  n <- nrow(D)
  qr_resid <- switch(residuals, full = design$qr, restricted = design$qr0)

  # d: number of coefficients of the model used for the residuals
  d <- qr_resid$rank
  d_expected <- switch(residuals,
    full = ncol(design$modelmat),
    restricted = ncol(design$modelmat) - length(design$indexes_X))
  stopifnot(d == d_expected, n > d)

  E <- qr.resid(qr_resid, D)                   # n x m residuals
  if (small_sample_corr == "HC2") {
    lev <- switch(residuals, full = design$lev, restricted = design$lev0)
    E <- E * sqrt(.cit_hc2_weights(lev, n, d))
  }
  G <- t(design$H)                             # n x K, row i = gamma_i

  U <- do.call(cbind, lapply(seq_len(ncol(G)),
    FUN = function(k) {
      G[, k] * E
    }))

  if (ncol(U) <= n) {
    S <- crossprod(U)
  } else {
    S <- tcrossprod(U)
  }

  ev <- eigen(S / n, symmetric = TRUE, only.values = TRUE)$values

  if (small_sample_corr == "HC1") {
    ev <- ev * n / (n - d)
  }

  # Sigma_hat is Positive Semi-Definite: clip round-off negatives and drop the
  # numerically null part of the spectrum (does not contribute to the chi-square
  # mixture).
  return(ev[ev > max(ev) * 1e-10])

}

# Eigenvalues of the factorized covariance (H H'/n) %x% Cov(D) used until
# citcdf 1.1.x (variance = "independent"). Only valid when the
# indicators D are independent of both X and Z. Within a gene it equals the former closed form min(p_j, p_j') - p_j p_j', and it also supplies the between-gene blocks.
.cit_independent_ev <- function(D, design) {
  n <- nrow(D)
  Dc <- D - matrix(colMeans(D), nrow = n, ncol = ncol(D), byrow = TRUE)
  ev_H <- eigen(tcrossprod(design$H) / n, symmetric = TRUE, only.values = TRUE)$values
  ev_D <- eigen(crossprod(Dc) / n, symmetric = TRUE, only.values = TRUE)$values

  return(as.vector(outer(ev_H, ev_D)))
}

# Warns when a level of a factor in X has fewer than 10 observations.
.cit_check_small_levels <- function(X) {
  X <- as.data.frame(X)
  counts <- unlist(lapply(X, FUN = function(x) {
    if (is.factor(x) || is.character(x)) {
      return(min(table(droplevels(as.factor(x)))))
    }
    return(NULL)
  }))

  if (length(counts) > 0 && min(counts) < 10) {
    warning("The asymptotic test is anti-conservative when a level of X has ",
      "few observations (here ", min(counts), "), regardless of n. Consider ",
      "residuals = \"restricted\" (conservative).", call. = FALSE)
  }

  return(invisible(NULL))
}

# Warns when the sandwich estimator is unreliable: n small relative to the
# number d of coefficients of the model used for the residuals (d/n > 0.05 for
# a single outcome with full residuals, d/n > 0.1 otherwise), or single
# outcome with full residuals and n < 30. The permutation test is only
# suggested with few coefficients. Messages when a correction is used for a
# gene set without need.
.cit_check_small_sample <- function(n, d, small_sample_corr, gene_set = FALSE,
                                    restricted = FALSE) {
  corrected <- small_sample_corr != "none"
  needed <- d / n > ifelse(gene_set || restricted, 0.1, 0.05)
  tiny <- n < 30 && !gene_set && !restricted
  if (d <= 3) {
    perm <- " Consider the permutation test."
  } else {
    perm <- ""
  }

  if (!corrected && (needed || tiny)) {
    warning("The asymptotic test is anti-conservative with n = ", n,
      " observations for ", d, " model coefficients. Consider ",
      "small_sample_corr = \"HC2\".", call. = FALSE)
  } else if (corrected && tiny) {
    warning("The asymptotic test can remain anti-conservative with n = ", n,
      " < 30, even with small_sample_corr = \"", small_sample_corr, "\".",
      perm, call. = FALSE)
  } else if (corrected && gene_set && !needed) {
    message("small_sample_corr = \"", small_sample_corr, "\" is conservative ",
      "for gene sets and reduces power: it is only needed when d/n > 0.1 ",
      "(here d = ", d, ", n = ", n, ").")
  }

  return(invisible(NULL))
}

# HC2 weights 1 / (1 - leverage), replaced by n / (n - d) for a leverage of 1.
.cit_hc2_weights <- function(lev, n, d) {
  w <- rep(n / (n - d), length(lev))
  below_1 <- 1 - lev > sqrt(.Machine$double.eps)
  w[below_1] <- 1 / (1 - lev[below_1])
  return(w)
}

# Warns when observations have a leverage of 1.
.cit_check_leverage <- function(lev) {
  n_lev1 <- sum(1 - lev <= sqrt(.Machine$double.eps))
  if (n_lev1 > 0) {
    warning(n_lev1, " observation(s) with a leverage of 1 (e.g. alone in a ",
      "level of a factor): their HC2 weight is replaced by n/(n-d), and the ",
      "asymptotic test is anti-conservative.", call. = FALSE)
  }

  return(invisible(NULL))
}

# Checks of the sandwich estimator, run once per call.
.cit_check_sandwich <- function(n, design, X, residuals, small_sample_corr,
                                gene_set = FALSE) {
  if (residuals == "full") {
    .cit_check_small_sample(n, design$qr$rank, small_sample_corr,
      gene_set = gene_set)
    .cit_check_small_levels(X)
    lev <- design$lev
  } else {
    .cit_check_small_sample(n, design$qr0$rank, small_sample_corr,
      gene_set = gene_set, restricted = TRUE)
    lev <- design$lev0
  }
  if (small_sample_corr == "HC2") {
    .cit_check_leverage(lev)
  }

  return(invisible(NULL))
}

.cit_warn_independent <- function() {
  warning("variance = \"independent\" is deprecated and will be removed in a ",
    "future release.\n  It is only valid when the outcome is independent of ",
    "both X and Z, and is very conservative when Z affects the outcome. ",
    "Use variance = \"sandwich\" (the default).", call. = FALSE)
}

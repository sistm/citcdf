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
  from <- if (length(which(Y == 0)) == 0) min(Y) else min(Y[-which(Y == 0)])
  seq(from = from, to = max(Y), length.out = number_y)
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
  list(modelmat = modelmat, indexes_X = indexes_X, H = H, XtXinv_X = XtXinv_X)
}

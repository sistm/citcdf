#' Shared internal utility functions
#'
#' Ddefinition of the y-threshold grid used by ccdf(), cit_asymp(), cit_perm()
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

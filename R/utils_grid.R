# Shared definition of the y-threshold grid.
#
# Used by cit_asymp(), cit_perm() and cit_gsa(), so that the
# asymptotic and permutation tests are evaluate the process on
# the same thresholds. Callers always drop the last grid point (`y[-p]`), so
# `to = max(Y)` is the first threshold NOT used: the degenerate threshold at
# which the indicator is identically 1 is already excluded by the loop.
#
# @keywords internal
# @noRd
.cit_y_grid <- function(Y, space_y, number_y) {
  Y <- as.numeric(Y)
  if (!space_y) {
    return(sort(unique(Y)))
  }
  from <- if (length(which(Y == 0)) == 0) min(Y) else min(Y[-which(Y == 0)])
  seq(from = from, to = max(Y), length.out = number_y)
}

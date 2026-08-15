#' Compute the conditional permutations
#'
#' @param X a numeric or factor vector of length \code{n}
#' containing the variable to be tested (the condition to be tested).
#' Multi-variables \code{X} are supported if \code{Z} is \code{NULL}.
#'
#' @param Z a numeric vector of length \code{n}
#' containing the covariate. Multiple variables are not allowed.
#'
#' @param n_perm the number of permutations. Default is \code{100}.
#'
#' @return a list with the permuted label vector for each permutation
#'
#' @details The permutations are conditional on \code{Z}:
#' \itemize{
#'  \item When \code{Z} is \code{NULL}, whole rows of \code{X} are permuted
#'  uniformly.
#'  \item When \code{Z} is a factor, integer, logical or character, \code{X} is
#'  permuted independently within each stratum of \code{Z}: the empirical
#'  association between \code{X} and \code{Z} is preserved exactly.
#'  \item When \code{Z} is continuous, no exact strata exist and the draw is
#'  delegated to \code{\link{perm_cont}}, which matches each observation to a
#'  swap of similar fitted value of \code{X} given \code{Z}.
#' }
#'
#' @seealso \code{\link{perm_cont}}, \code{\link{cit_perm}}
#'
#' @export
#'
#' @examples
#' set.seed(123)
#' X <- rbinom(n = 100, size = 1, prob = 0.5)
#' Z <- rnorm(100, 0, 1)
#' X_perm(data.frame(X), data.frame(Z), 100)
#'
X_perm <- function(X, Z, n_perm = 100) {

  stopifnot(is.data.frame(X))
  stopifnot(is.data.frame(Z) | is.null(Z))
  if (!is.null(Z)) {
    stopifnot(ncol(X) < 2)
  }
  stopifnot(ncol(Z) < 2 | is.null(Z))


  # pas besoin pas de Y
  # stopifnot(nrow(X) == n)
  # stopifnot(nrow(Z) == n | is.null(Z))


  if (is.null(Z)) {
    colnames(X) <- paste0("X", seq_len(ncol(X)))
  } else { # with covariates Z
    colnames(X) <- paste0("X", seq_len(ncol(X)))
    colnames(Z) <- paste0("Z", seq_len(ncol(Z)))
  }

  X_star <- list()

  if (is.null(Z)) {
    for (k in seq_len(n_perm)) {
      X_star[[k]] <- X[sample(nrow(X)), , drop = FALSE]
      rownames(X_star[[k]]) <- NULL
    }

  } else {

    sample_X <- function(z) {
      idx <- seq_along(z)
      for (zj in unique(z)) {
        s <- which(z == zj)
        idx[s] <- s[sample.int(length(s))]
      }
      return(idx)
    }

    xcol <- X[, 1]
    x <- as.numeric(xcol) # perm_cont() needs a numeric X

    for (k in seq_len(n_perm)) {
      z <- Z[, 1]
      if (is.factor(z) || is.integer(z) || is.logical(z) || is.character(z)) {
        # discrete Z: permute X within each stratum of Z
        zz <- as.numeric(as.factor(z))
        xs <- xcol[sample_X(zz)]
      } else if (is.numeric(z)) {
        # continuous Z: distance-weighted conditional permutation
        codes <- perm_cont(X = x, Z = z)
        if (is.factor(xcol)) {
          xs <- factor(levels(xcol)[codes], levels = levels(xcol))
        } else {
          xs <- codes
        }
      } else {
        stop("unsupported class for 'Z': ", paste(class(z), collapse = "/"),
          ". 'Z' must be a factor, integer, logical, character or numeric column.")
      }
      X_star[[k]] <- if (is.factor(xcol)) data.frame(X = xs) else xs
    }
  }

  return(X_star)
}

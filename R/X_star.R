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
#' @seealso \code{\link{perm_cont}}, \code{\link{cit_perm}}
#'
#' @export
#'
#' @examples
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

    sample_X <- function(X, Z, z) {
      X_sampled <- rep(NA, length(Z))
      for (zj in unique(z)) {
        X_sampled[Z == zj] <- sample(X[Z == zj])
      }
      return(X_sampled)
    }

    for (k in seq_len(n_perm)) {
      z <- Z[, 1]
      x <- if (is.factor(X[, 1])) as.numeric(levels(X[, 1]))[X[, 1]] else as.numeric(X[, 1])
      X_star[[k]] <- if (is.factor(z) || is.integer(z) || is.logical(z) ||
        is.character(z)) {
        # discrete Z: permute X within each stratum of Z
        zz <- as.numeric(as.factor(z))
        sample_X(X[, 1], zz, unique(zz))
      } else if (is.numeric(z)) {
        # continuous Z: distance-weighted conditional permutation
        perm_cont(X = x, Z = z)
      } else {
        stop("unsupported class for 'Z': ", paste(class(z), collapse = "/"),
          ". 'Z' must be a factor, integer, logical, character or numeric column.")
      }
      if (is.factor(X[, 1])) {
        X_star[[k]] <- data.frame(X = as.factor(X_star[[k]]))
      }

      # colnames(X_star) <- sapply(seq_len(ncol(X)), function(i){paste0('X',i)})
    }
  }

  return(X_star)
}

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
      X_star[[k]] <- switch(class(Z[, 1]),
        "factor" = sample_X(X[, 1], as.numeric(Z[, 1]),
          unique(as.numeric(Z[, 1]))), # to be double-checked
        "integer" = sample_X(X[, 1], as.numeric(Z[, 1]),
          unique(as.numeric(Z[, 1]))), # to be double-checked
        "numeric" = perm_cont(X = if (is.factor(X[, 1])) {
          as.numeric(levels(X[, 1]))[X[, 1]]
        } else {
          as.numeric(X[, 1])
        }, Z = Z[, 1])
      )
      if (is.factor(X[, 1])) {
        X_star[[k]] <- data.frame(X = as.factor(X_star[[k]]))
      }

      # colnames(X_star) <- sapply(seq_len(ncol(X)), function(i){paste0('X',i)})
    }
  }

  return(X_star)
}

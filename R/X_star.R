#' Compute the conditionnal permutation
#'
#' @param X a numeric or factor vector of length \code{n}
#' containing the variable to be tested (the condition to be tested).
#'
#' @param Z a numeric vector of length \code{n}
#' containing the covariate. Multiple variables are not allowed.
#'
#' @param n_perm the number of permutations. Default is \code{100}.
#'
#' @returns a list with the permuted label vector for each permutation
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
  stopifnot(ncol(X) < 2)
  stopifnot(ncol(Z) < 2 | is.null(Z))


  # pas besoin pas de Y
  # stopifnot(nrow(X) == n)
  # stopifnot(nrow(Z) == n | is.null(Z))


  if (is.null(Z)) {
    colnames(X) <- sapply(1:ncol(X), function(i) {
      paste0("X", i)
    })
  } else { # with covariates Z
    colnames(X) <- sapply(1:ncol(X), function(i) {
      paste0("X", i)
    })
    colnames(Z) <- sapply(1:ncol(Z), function(i) {
      paste0("Z", i)
    })
  }

  X_star <- list()

  if (is.null(Z)) {
    for (k in 1:n_perm) {
      X_star[[k]] <- data.frame(X = X[sample(1:nrow(X)), ])
      # colnames(X_star) <- sapply(1:ncol(X), function(i){paste0('X',i)})
    }

  } else {

    sample_X <- function(X, Z, z) {
      X_sampled <- rep(NA, length(Z))
      for (zj in unique(z)) {
        X_sampled[Z == zj] <- sample(X[Z == zj])
      }
      return(X_sampled)
    }

    for (k in 1:n_perm) {
      X_star[[k]] <- switch(class(Z[, 1]),
        "factor" = sample_X(X[, 1], as.numeric(Z[, 1]), unique(as.numeric(Z[, 1]))), # vérifier
        "integer" = sample_X(X[, 1], as.numeric(Z[, 1]), unique(as.numeric(Z[, 1]))), # vérifier
        "numeric" = perm_cont(X = if (is.factor(X[, 1])) {
          as.numeric(levels(X[, 1]))[X[, 1]]
        } else {
          as.numeric(X[, 1])
        }, Z = Z[, 1])
      )
      if (is.factor(X[, 1])) {
        X_star[[k]] <- data.frame(X = as.factor(X_star[[k]]))
      }

      # colnames(X_star) <- sapply(1:ncol(X), function(i){paste0('X',i)})
    }
  }

  return(X_star)
}

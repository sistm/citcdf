#' Permutation test for conditional independence
#'
#' @param Y a numeric vector of length \code{n} to test for conditional independence
#' with \code{X} adjusted on \code{Z}
#'
#' @param X a data frame of size \code{n x p} of numeric or factor vector(s)
#' containing the variable(s) to be tested for conditional independence
#' against \code{X} adjusted on \code{Z}. Multi-variables \code{X} are
#' supported if \code{Z} is \code{NULL}.
#'
#' @param Z a data.frame of size \code{n x 1} of numeric or factor vector
#' containing the covariate to condition the independence
#' test upon. Multiple covariates are not supported for permutation.
#'
#' @param X_star a list of \code{n_perm} permuted designs, as returned
#' by \code{\link{X_perm}}. Default is \code{NULL}, in which case
#' \code{X_perm(X, Z, n_perm = n_perm)} is called internally. Supply it
#' explicitly when several outcomes must be scored against the same
#' permutations, or to avoid redrawing them inside a loop; see \emph{Details}.
#'
#' @param n_perm the number of permutations. Default is \code{100}. When
#' \code{X_star} is supplied it must hold at least \code{n_perm} elements;
#' only the first \code{n_perm} are used.
#'
#' @param space_y a logical flag indicating whether the y thresholds are spaced out.
#' When \code{space_y} is \code{TRUE}, a regular sequence between the minimum and
#' the maximum of the observations is used. Default is \code{FALSE}.
#'
#' @param number_y an integer value indicating the number of y thresholds (and therefore
#' the number of regressions) to perform the test. Only used if \code{space_y}
#' is \code{TRUE}. Default is \code{10}.
#'
#' @details The \code{space_y} / \code{number_y} grid controls both the
#' resolution of the statistic and its computational cost. See
#' \code{\link{cit_multi}} for details on this trade-off.
#'
#' Leaving \code{X_star} as \code{NULL} is the convenient form for a single
#' outcome. Across several outcomes it is not equivalent to supplying one:
#' each call would draw its own permutations, whereas \code{\link{cit_multi}}
#' and \code{\link{cit_gsa}} deliberately build one pool with
#' \code{\link{X_perm}} and reuse it for every gene, so that all genes are
#' scored against the same permuted designs. Pass a shared \code{X_star} if
#' you are looping over outcomes yourself.
#'
#' @seealso \code{\link{perm_cont}}, \code{\link{X_perm}}, \code{\link{cit_multi}}
#'
#' @export
#'
#' @return A data frame with the following elements:
#' \itemize{
#'   \item \code{score} contains the number of permutations whose test
#'   statistic is greater than or equal to the observed one.
#'   \item \code{raw_pval} contains the raw p-values for a given gene computed
#'   from \code{n_perm} permutations.
#'   \item \code{test_statistic} contains the observed test statistic for a given
#'   gene. It is the same quantity returned by \code{\link{cit_asymp}}.
#' }
#'
#' @references Gauthier M, Agniel D, Thiébaut R & Hejblum BP (2021).
#' Distribution-free complex hypothesis testing for single-cell RNA-seq
#' differential expression analysis, \emph{bioRxiv} 445165.
#' \doi{10.1101/2021.05.21.445165}.
#'
#' @examples
#'
#' set.seed(123)
#' X <- data.frame(X = as.factor(rbinom(n = 100, size = 1, prob = 0.5)))
#' Y <- (X$X == 1) * rnorm(100) + (X$X == 0) * rnorm(100, mean = 0.5)
#'
#' # the permuted designs are drawn internally when X_star is left NULL
#' res_perm <- cit_perm(Y, X, n_perm = 10)
#' res_perm
#'
#' # supplying them explicitly is equivalent, and is what to do when several
#' # outcomes must be scored against the same permutations
#' X_star <- X_perm(X, Z = NULL, n_perm = 10)
#' res_perm_shared <- cit_perm(Y, X, X_star = X_star, n_perm = 10)
#'
#' # adjusting for a covariate Z
#' Z <- data.frame(Z = rnorm(100))
#' res_perm_adj <- cit_perm(Y, X, Z = Z, n_perm = 10)
#' res_perm_adj
cit_perm <- function(Y, X, Z = NULL, X_star = NULL, n_perm = 100, space_y = FALSE, number_y = 10) {

  stopifnot(is.vector(Y))
  stopifnot(is.data.frame(X))
  stopifnot(is.data.frame(Z) | is.null(Z))
  if (!is.null(Z)) {
    stopifnot(ncol(X) < 2)
  }
  stopifnot(ncol(Z) < 2 | is.null(Z))

  n <- length(Y)
  stopifnot(nrow(X) == n)
  stopifnot(nrow(Z) == n | is.null(Z))
  .cit_check_Y(Y)


  if (is.null(X_star)) {
    X_star <- X_perm(X, Z, n_perm = n_perm)
  }
  stopifnot(is.list(X_star))
  if (length(X_star) < n_perm) {
    stop("'X_star' holds ", length(X_star), " permuted design(s) but 'n_perm' ",
      "is ", n_perm, ". Supply at least 'n_perm' of them, or leave 'X_star' ",
      "as NULL to have them drawn here.", call. = FALSE)
  } else if ((length(X_star) > n_perm)) {
    warning("'X_star' holds ", length(X_star), " permuted design(s), which is
            larger than 'n_perm', which itself is ", n_perm, ". Only the first ",
      n_perm, "values from 'X_star' are used.", call. = FALSE)
  }



  if (is.null(Z)) {
    colnames(X) <- paste0("X", seq_len(ncol(X)))
    modelmat <- model.matrix(~., data = X)
  } else { # with covariates Z
    colnames(X) <- paste0("X", seq_len(ncol(X)))
    colnames(Z) <- paste0("Z", seq_len(ncol(Z)))
    modelmat <- model.matrix(~., data = cbind(X, Z))
  }

  indexes_X <- which(substring(colnames(modelmat), 1, 1) == "X")

  n_Y_all <- length(Y)
  H <- n_Y_all * (solve(crossprod(modelmat)) %*% t(modelmat))[indexes_X, , drop = FALSE]

  Y <- as.numeric(Y)
  oY <- order(Y)
  y <- .cit_y_grid(Y, space_y, number_y)
  p <- length(y) # number of thresholds used

  index_jumps <- findInterval(y[-p], Y[oY])
  beta <- c(apply(X = H[, oY, drop = FALSE], MARGIN = 1, FUN = cumsum)[index_jumps, , drop = FALSE]) / n_Y_all
  test_stat_obs <- sum(beta^2) * n_Y_all


  test_stat_perm <- rep(NA, n_perm)

  if (is.null(Z)) {
    XtXinv_X <- solve(crossprod(modelmat))[indexes_X, , drop = FALSE]
    for (k in seq_len(n_perm)) {
      modelmat_perm <- model.matrix(~., data = X_star[[k]])

      H_perm <- n_Y_all * tcrossprod(XtXinv_X, modelmat_perm)
      beta_perm <- c(apply(X = H_perm[, oY, drop = FALSE], MARGIN = 1, FUN = cumsum)[index_jumps, , drop = FALSE]) / n_Y_all

      test_stat_perm[k] <- sum(beta_perm^2) * n_Y_all
    }
  } else {

    for (k in seq_len(n_perm)) {
      # colnames(X_star) <- sapply(seq_len(ncol(X)), function(i){paste0('X',i)})
      modelmat_perm <- model.matrix(~., data = cbind(X_star[[k]], Z))

      H_perm <- n_Y_all * (solve(crossprod(modelmat_perm)) %*% t(modelmat_perm))[indexes_X, , drop = FALSE]
      beta_perm <- c(apply(X = H_perm[, oY, drop = FALSE], MARGIN = 1, FUN = cumsum)[index_jumps, , drop = FALSE]) / n_Y_all

      test_stat_perm[k] <- sum(beta_perm^2) * n_Y_all
    }
  }

  score <- sum(test_stat_perm >= test_stat_obs)
  pval <- (score + 1) / (n_perm + 1)
  return(data.frame(score = score, raw_pval = pval,
    test_statistic = test_stat_obs))

}

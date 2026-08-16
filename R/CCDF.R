#' Estimate the empirical conditional cumulative distribution function
#'
#' @param Y a numeric vector of size \code{n} containing the
#' preprocessed expressions from \code{n} samples (or cells).
#'
#' @param X a data frame containing numeric or factor vector(s) of size \code{n}
#' containing the variable(s) to be tested (the condition(s) to be tested).
#'
#' @param Z a data frame containing numeric or factor vector(s) of size \code{n}
#' containing the covariate(s).
#'
#' @param method a character string indicating which method to use to
#' compute the CCDF, either \code{'OLS'} or \code{'logistic'}.
#' Default is \code{'OLS'} for greater computational speed.
#'
#' @param fast a logical flag indicating whether the fast implementation of
#' logistic regression should be used. Only if \code{method == 'logistic'}.
#' Default is \code{TRUE}.
#'
#' @param space_y a logical flag indicating whether the y thresholds are spaced.
#' When \code{space_y} is \code{TRUE}, a regular sequence between the minimum and
#' the maximum of the observations is used. Default is \code{FALSE}.
#'
#' @param number_y an integer value indicating the number of y thresholds (and therefore
#' the number of regressions) to perform the test. Default is \code{length(Y)}.
#'
#' @importFrom stats model.matrix
#'
#' @details The \code{space_y} / \code{number_y} grid controls both the
#' resolution of the statistic and its computational cost. See
#' \code{\link{cit_multi}} for details on this trade-off.
#'
#' @references Gauthier M, Agniel D, Thiébaut R & Hejblum BP (2021).
#' Distribution-free complex hypothesis testing for single-cell RNA-seq
#' differential expression analysis, *bioRxiv* 445165.
#' [DOI: 10.1101/2021.05.21.445165](https://doi.org/10.1101/2021.05.21.445165).
#'
#' @seealso \code{\link{plot_compare_ccdf}}, \code{\link{cit_asymp}}
#'
#' @export
#'

#' @return A list with the following elements:\itemize{
#'   \item \code{cdf}: a vector of the cumulative distribution function of a given gene.
#'   \item \code{ccdf}: a vector of the conditional cumulative distribution function of a given gene, computed
#'   given \code{X}. Only if \code{Z} is \code{NULL}.
#'   \item \code{ccdf_nox}: a vector of the conditional cumulative distribution function of a given gene, computed
#'   given \code{Z} only (i.e. \code{X} is ignored.). Only if \code{Z} is not \code{NULL}.
#'   \item \code{ccdf_x}: a vector of the conditional cumulative distribution function of a given gene, computed
#'   given \code{X} and \code{Z}. Only if \code{Z} is not \code{NULL}.
#'   \item \code{y_sort}: a vector of the sorted expression points at which the CDF and the CCDFs are calculated.
#'   \item \code{x_sort}: a vector of the variables associated with \code{y_sort}.
#'   \item \code{z_sort}: a vector of the covariates associated with \code{y_sort}. Only if \code{Z} is not \code{NULL}.
#' }
#'
#' @examples
#'
#' set.seed(123)
#' n <- 500
#' X <- as.factor(rbinom(n = n, size = 1, prob = 0.5))
#' Y <- ((X == 1) * rnorm(n = n, 0, 1)) + ((X == 0) * rnorm(n = n, 0.5, 1))
#' res <- ccdf(Y, data.frame(X = X), method = "OLS")
ccdf <- function(Y, X, Z = NULL, method = c("OLS", "logistic"),
                 fast = TRUE, space_y = FALSE, number_y = length(Y)) {

  if (length(method) > 1) {
    method <- method[1]
  }
  stopifnot(method %in% c("OLS", "logistic"))

  # --- structural checks on the inputs ----------------------------------------
  # Y[i] is the response for row i of X (and of Z), so the three must line up.
  # model.matrix() silently drops rows containing NA, which would shift X out of
  # step with Y and produce a wrong answer rather than an error, so NAs in the
  # design are refused up front.
  if (!is.data.frame(X)) {
    stop("'X' must be a data frame with one column per variable to condition on, ",
      "not a ", class(X)[1], ". Did you mean data.frame(X = X)?")
  }
  if (!is.null(Z) && !is.data.frame(Z)) {
    stop("'Z' must be a data frame or NULL, not a ", class(Z)[1],
      ". Did you mean data.frame(Z = Z)?")
  }

  n_obs <- length(Y)
  if (nrow(X) != n_obs) {
    stop("'X' has ", nrow(X), " row(s) but 'Y' has ", n_obs, " observation(s). ",
      "'Y[i]' is the response for row i of 'X', so the two must have the same length.")
  }
  if (!is.null(Z) && nrow(Z) != n_obs) {
    stop("'Z' has ", nrow(Z), " row(s) but 'Y' has ", n_obs, " observation(s). ",
      "'Y[i]' is the response for row i of 'Z', so the two must have the same length.")
  }
  if (anyNA(X)) {
    stop("'X' contains ", sum(is.na(X)), " NA value(s). ",
      "Rows with NA are silently dropped when the model matrix is built, which ",
      "would misalign 'X' with 'Y'. Please remove or impute them before calling ccdf().")
  }
  if (!is.null(Z) && anyNA(Z)) {
    stop("'Z' contains ", sum(is.na(Z)), " NA value(s). ",
      "Rows with NA are silently dropped when the model matrix is built, which ",
      "would misalign 'Z' with 'Y'. Please remove or impute them before calling ccdf().")
  }


  if (!is.numeric(Y)) {
    warning("Converting Y to a numeric vector.\n",
      "This should have been done beforehand.")
    Y <- as.numeric(Y)
  }

  if (anyNA(Y)) {
    keep <- !is.na(Y)
    warning("`Y` contains ", sum(!keep), " NA values. ",
      "\nCurrently those are ignored in the computations ",
      "(dropped from `Y`, `X` and `Z`), but ",
      "you should think carefully about where do those NA/NaN ",
      "come from...")
    Y <- Y[keep]
    X <- X[keep, , drop = FALSE]
    if (!is.null(Z)) Z <- Z[keep, , drop = FALSE]
  }

  y <- .cit_y_grid(Y, space_y, number_y)

  output <- NULL

  if (is.null(Z)) {
    n_Y <- length(Y)
    # temp_order <- sort(Y,index.return=TRUE)$ix
    # y <- sort(unique(Y))
    # y_sort <- sort(Y)
    # x_sort <- X[temp_order]
    # modelmat <- model.matrix(Y~X)


    colnames(X) <- paste0("X", seq_len(ncol(X)))
    modelmat <- model.matrix(~., data = X)


    ind_X <- which(substring(colnames(modelmat), 1, 1) == "X")
    P <- solve(crossprod(modelmat)) %*% t(modelmat)

    cdf <- list()
    ccdf <- list()
    x_sort <- list()
    y_sort <- list()

    for (i in seq_len(length(y) - 1)) {
      Ylow <- Y <= y[i]
      if (i == 1) {
        w <- Ylow
      } else {
        w <- Ylow & (Y > y[i - 1])
      }
      x_sort[[i]] <- X[w, ]
      y_sort[[i]] <- Y[w]
      indi_Y <- 1 * Ylow

      # unCDF
      cdf[[i]] <- rep(sum(indi_Y) / n_Y, sum(w))

      if (length(unique(indi_Y)) == 1) {
        ccdf[[i]] <- rep(1, sum(w))
      } else {

        if (method == "logistic") {
          # CDF
          if (fast) {
            # fast
            fit <- RcppNumerical::fastLR(x = modelmat, y = indi_Y,
              eps_f = 1e-08, eps_g = 1e-08)
          } else {
            fit <- glm(indi_Y ~ ., data = cbind.data.frame(indi_Y, X),
              family = binomial(link = "logit"))
          }
          # fastLR() and glm() both return the fitted probabilities on the full design
          ccdf[[i]] <- fit$fitted.values[w]
        } else if (method == "OLS") {
          coefs_ols <- P %*% indi_Y
          ccdf[[i]] <- (modelmat[w, , drop = FALSE] %*% coefs_ols)[, 1]
        }
      }
    }


    ccdf_unlisted <- unlist(ccdf, use.names = FALSE)
    cdf_unlisted <- unlist(cdf, use.names = FALSE)
    output <- list(cdf = cdf_unlisted, ccdf = ccdf_unlisted, y = unlist(y_sort), x = unlist(x_sort))
    class(output) <- "ccdf"

  } else {
    n_Y <- length(Y)

    colnames(X) <- paste0("X", seq_len(ncol(X)))
    colnames(Z) <- paste0("Z", seq_len(ncol(Z)))
    modelmat <- model.matrix(~., data = cbind.data.frame(X, Z))

    ind_X <- which(substring(colnames(modelmat), 1, 1) == "X")
    mm_nox <- modelmat[, -ind_X, drop = FALSE]
    P <- solve(crossprod(modelmat)) %*% t(modelmat)
    P_nox <- solve(crossprod(mm_nox)) %*% t(mm_nox)

    x_sort <- list()
    y_sort <- list()
    z_sort <- list()

    ccdf_x <- list()
    ccdf_nox <- list()
    cdf <- list()

    for (i in seq_len(length(y) - 1)) {
      # new_data <- data.frame(X[w],Z[w])
      # names(new_data) <- c("X","Z")

      # cdf[[i]] <- rep(sum(indi_Y)/n_Y, sum(w))

      Ylow <- Y <= y[i]

      if (i == 1) {
        w <- Ylow
      } else {
        w <- Ylow & (Y > y[i - 1])
      }

      x_sort[[i]] <- X[w, ]
      y_sort[[i]] <- Y[w]
      z_sort[[i]] <- Z[w, ]
      indi_Y <- 1 * Ylow

      # unCDF
      cdf[[i]] <- rep(sum(indi_Y) / n_Y, sum(w))

      if (length(unique(indi_Y)) == 1) {
        ccdf_x[[i]] <- rep(indi_Y[1], sum(w))
        ccdf_nox[[i]] <- rep(indi_Y[1], sum(w))

      } else if (method == "logistic") {
        mm_nox <- modelmat[, -ind_X, drop = FALSE]
        if (fast) {
          fit_x   <- RcppNumerical::fastLR(x = modelmat, y = indi_Y)
          fit_nox <- RcppNumerical::fastLR(x = mm_nox,   y = indi_Y)
        } else {
          fit_x   <- glm.fit(x = modelmat, y = indi_Y, family = binomial())
          fit_nox <- glm.fit(x = mm_nox,   y = indi_Y, family = binomial())
        }
        ccdf_x[[i]]   <- fit_x$fitted.values[w]
        ccdf_nox[[i]] <- fit_nox$fitted.values[w]

      } else if (method == "OLS") {
        coefs_ols_x <- P %*% indi_Y
        ccdf_x[[i]] <- (modelmat[w, , drop = FALSE] %*% coefs_ols_x)[, 1]
        coefs_ols_nox <- P_nox %*% indi_Y
        ccdf_nox[[i]] <- (mm_nox[w, , drop = FALSE] %*% coefs_ols_nox)[, 1]
      }
    }

    ccdf_x_unlisted <- unlist(ccdf_x, use.names = FALSE)
    ccdf_nox_unlisted <- unlist(ccdf_nox, use.names = FALSE)
    cdf_unlisted <- unlist(cdf, use.names = FALSE)

    output <-  list(cdf = cdf_unlisted, ccdf_nox = ccdf_nox_unlisted, ccdf_x = ccdf_x_unlisted, y = unlist(y_sort), x = unlist(x_sort), z = unlist(z_sort))
    class(output) <- "ccdf"
  }


  return(output)
}

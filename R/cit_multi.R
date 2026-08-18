#' Multiple conditional independence testing
#'
#' @param M a \code{data.frame} or a \code{matrix} of size \code{n x r}
#' containing the different Y variables to test for conditional independence
#' with \code{X} adjusted on \code{Z}
#'
#' @param X a data frame of size \code{n x p} of numeric or factor vector(s)
#' containing the variable(s) to be tested for conditional independence
#' against \code{X} adjusted on \code{Z}. Multiple variables (\code{p>1})
#' are supported by the asymptotic test, and also by the permutation when
#' \code{Z} is \code{NULL}.
#'
#' @param Z a data frame of size \code{n x q} of numeric or factor vector(s)
#' containing the covariate(s) to condition the independence
#' test upon. Multiple covariates (\code{q>1}) are only supported by the
#' asymptotic test.
#'
#' @param test a character string indicating whether the \code{'asymptotic'} or
#' the \code{'permutation'} test is computed.
#' Default is \code{'asymptotic'}.
#'
#' @param n_perm the number of permutations. Default is \code{100}. Only used if
#' \code{test == 'permutation'}.
#'
#' @param adaptive a logical flag indicating whether adaptive additional permutations
#' should be performed. Default is \code{TRUE}. Only used if
#' \code{test == 'permutation'}.
#' Note \code{cit_gsa()} defaults to \code{adaptive = FALSE} whereas
#' \code{cit_multi()} defaults to \code{TRUE}.
#'
#' @param n_perm_adaptive a vector of the increasing numbers of
#' adaptive permutations to be performed when \code{adaptive} is \code{TRUE}
#' if p-values are below \code{thresholds}.
#' \code{length(n_perm_adaptive)} should be equal to \code{length(thresholds)+1}.
#' Default is \code{c(n_perm, n_perm, n_perm*3, n_perm*5)}.
#'
#' @param thresholds a vector of the decreasing thresholds to compute
#' adaptive permutations when \code{adaptive} is \code{TRUE}.
#' \code{length(thresholds)} should be equal to \code{length(n_perm_adaptive)-1}.
#' Default is \code{c(0.1, 0.05, 0.01)}.
#'
#'
#' @param parallel a logical flag indicating whether parallel computation
#' should be enabled. Default is \code{TRUE} if \code{interactive()} is
#' \code{TRUE}, else is \code{FALSE}.
#'
#' @param n_cpus an integer indicating the number of cores to be used for the
#' computations. Default is
#' \code{max(1L, paralleldetectCores(logical = FALSE) - 1L, na.rm = TRUE)}.
#' If \code{n_cpus = 1}, then sequential computations are used without any
#' parallelization.
#'
#' @param space_y a logical flag indicating whether the y thresholds are spaced out.
#' When \code{space_y} is \code{TRUE}, a regular sequence between the minimum and
#' the maximum of the observations is used. If \code{FALSE}, each unique
#' observed expression value is used as a distinct threshold. Default is \code{TRUE}.
#'
#' @param number_y an integer value indicating the number of y thresholds (and therefore
#' the number of regressions) to perform the test. Only used if \code{space_y}
#' is \code{TRUE}. Default is \code{10}.
#'
#'
#' @import pbapply
#' @importFrom parallel makeCluster stopCluster detectCores
#'
#'
#' @return A list with the following elements:\itemize{
#'   \item \code{which_test}: a character string carrying forward the value of
#'   the '\code{test}' argument indicating which test was performed (either
#'   'asymptotic' or 'permutation').
#'   \item \code{n_perm}: an integer carrying forward the value of the
#'   '\code{n_perm}' argument or '\code{n_perm_adaptive}' indicating the number of permutations performed
#'   (\code{NA} if asymptotic test was performed).
#'   \item \code{pvals}: computed p-values. A data frame with one row for
#'   each gene, and with 2 columns: the first one '\code{raw_pval}' contains
#'   the raw p-values, the second one '\code{adj_pval}' contains the FDR adjusted p-values
#'   using Benjamini-Hochberg correction. When '\code{test == "asymptotic"}', a
#'   third column '\code{test_statistic}' contains the gene-wise test statistics.
#' }
#'
#' @details With \code{space_y = FALSE} the test statistic uses every distinct
#' observed in \code{Y}), but its computation cost represents
#' \code{length(unique(Y))} regressions per gene.
#' With \code{space_y = TRUE}, it uses instead a regular grid of \code{number_y}
#' points, trading resolution for a computational cost that is independent of
#' \code{n}. The grid runs from the smallest non-zero observation (a mass of
#' exact zeros, as in count data, does not consume grid points) to
#' \code{max(Y)}. Raising \code{number_y} brings p-values
#' closer towards their \code{space_y = FALSE} values.
#'
#' @seealso \code{\link{cit_asymp}}, \code{\link{cit_perm}}, \code{\link{ccdf}}
#'
#' @references Gauthier M, Agniel D, Thiébaut R & Hejblum BP (2021).
#' Distribution-free complex hypothesis testing for single-cell RNA-seq
#' differential expression analysis, \emph{bioRxiv} 445165.
#' \doi{10.1101/2021.05.21.445165}.
#'
#' @export
#'
#' @examples
#'
#'
#' set.seed(123)
#' Z <- as.factor(rbinom(n = 100, size = 1, prob = 0.5))
#' X <- as.numeric(Z) - 1  + rnorm(n = 100, sd = 1)
#' r <- 500
#' Y <- replicate(r, as.numeric(Z) - 1)
#' Y <- (Y == 1) * rnorm(n = 100 * r, 0, 1) + (Y == 0) * rnorm(n = 100 * r, 0.5, 1)
#' res_asymp_unadj <- cit_multi(M = data.frame(Y = Y),
#'   X = data.frame(X = X),
#'   test = "asymptotic", parallel = FALSE)
#' mean(res_asymp_unadj$pvals$raw_pval < 0.05)
#' hist(res_asymp_unadj$pvals$raw_pval)
#'
#' res_asymp_adj <- cit_multi(M = data.frame(Y = Y),
#'   X = data.frame(X = X),
#'   Z = data.frame(Z = Z),
#'   test = "asymptotic", parallel = FALSE)
#' mean(res_asymp_adj$pvals$raw_pval < 0.05)
#' hist(res_asymp_adj$pvals$raw_pval)
#'
#' n <- 100
#' r <- 500
#' Z1 <- rbinom(n, size = 1, prob = 0.5)
#' Z2 <- rnorm(n) # rbinom(n, size=1, prob=0.5) + rnorm(n, sd=0.05)
#' X1 <- Z2 + rnorm(n, sd = 0.2)
#' X2 <- rnorm(n)
#' cor(X1, Z2)
#' Y <- replicate(r, Z2) + rnorm(n * r, 0, 3)
#' range(cor(Y, Z2))
#' range(cor(Y, X2))
#' res_asymp_unadj <- cit_multi(M = data.frame(Y = Y),
#'   X = data.frame(X1 = X1, X2 = X2),
#'   test = "asymptotic", parallel = FALSE)
#' mean(res_asymp_unadj$pvals$raw_pval < 0.05)
#' hist(res_asymp_unadj$pvals$raw_pval)
#'
#' res_asymp_adj <- cit_multi(M = data.frame(Y = Y),
#'   X = data.frame(X1 = X1, X2 = X2),
#'   Z = data.frame(Z1 = Z1, Z2 = Z2),
#'   test = "asymptotic", parallel = FALSE)
#' mean(res_asymp_adj$pvals$raw_pval < 0.05)
#' hist(res_asymp_adj$pvals$raw_pval)
#'
#' # permutation test, on a subset of the genes to keep the example short
#' res_perm_unadj <- cit_multi(M = data.frame(Y = Y[, 1:20]),
#'   X = data.frame(X1 = X1),
#'   test = "permutation", adaptive = FALSE, n_perm = 50,
#'   parallel = FALSE)
#' mean(res_perm_unadj$pvals$raw_pval < 0.05)
#'
#' \donttest{
#' # adaptive permutations spend extra stages only on the smallest p-values
#' res_perm_adj <- cit_multi(M = data.frame(Y[, 1:20]), # data.frame(Y),
#'   X = data.frame(X = X),
#'   Z = data.frame(Z = Z),
#'   test = "permutation", n_perm = 50, # 2000,
#'   parallel = FALSE)
#' mean(res_perm_adj$pvals$raw_pval < 0.05)
#' }
cit_multi <- function(M,
                      X,
                      Z = NULL,
                      test = c("asymptotic", "permutation"),
                      n_perm = 100,
                      n_perm_adaptive = c(n_perm, n_perm, n_perm * 3, n_perm * 5),
                      thresholds = c(0.1, 0.05, 0.01),
                      parallel = interactive(),
                      n_cpus = max(1L, detectCores(logical = FALSE) - 1L, na.rm = TRUE),
                      adaptive = TRUE,
                      space_y = TRUE,
                      number_y = 10) {
  # check
  if (is.matrix(M)) {
    M <- as.data.frame(M)
  }

  stopifnot(is.data.frame(M))
  stopifnot(is.data.frame(X))
  stopifnot(is.data.frame(Z) | is.null(Z))
  stopifnot(is.logical(parallel))
  stopifnot(is.logical(adaptive))
  stopifnot(is.numeric(n_perm))

  if (anyNA(M)) {
    warning("'M' contains", sum(is.na(M)), "NA values. ",
      "\nCurrently they are ignored in the computations but ",
      "you should think carefully about where do those NA/NaN ",
      "come from...")
    M <- M[, colSums(is.na(M)) == 0, drop = FALSE]
  }

  M_colnames <- colnames(M)
  r <- ncol(M)
  n <- nrow(M)
  stopifnot(nrow(X) == n)
  stopifnot(nrow(Z) == n | is.null(Z))
  .cit_check_M(M)

  if (length(test) > 1) {
    test <- test[1]
  }
  stopifnot(test %in% c("asymptotic", "permutation"))

  if (test == "permutation") {

    if (adaptive) {
      if ((length(n_perm_adaptive) != (length(thresholds) + 1))) {
        warning("length of thresholds + 1 must be equal to length of n_perm_adaptive. \n",
          "Consider using the default parameters.")
      }
    } else {
      N_possible_perms <- factorial(n)
      if (n_perm > N_possible_perms) {
        warning("The number of permutations requested 'n_perm' is ",
          n_perm, "which is larger than the total number of ",
          "existing permutations ", N_possible_perms,
          ". Try a lower number for 'n_perm' (currently ",
          "running with 'nperm=", N_possible_perms, "').")
        n_perm <- N_possible_perms
      }
    }
  }

  if (space_y) {
    if (is.null(number_y)) {
      warning("Missing argument", number_y, ". No spacing is used.")
      space_y <- FALSE
    }
  }

  # parallel

  if (parallel) {

    if (.Platform$OS.type == "unix") {
      par_clust <- n_cpus
    } else {
      par_clust <- parallel::makeCluster(n_cpus)
    }

  } else {
    par_clust <- 1L
  }



  # Test ----
  ## permutations ----
  if (test == "permutation") {
    # multi-column X is supported without Z, conditional permutation scheme
    # (Z present) is restricted to a single column.
    if (!is.null(Z) && ncol(X) > 1) {
      stop("the permutation test supports a single column in 'X' only when 'Z' ",
        "is provided. Use test = \"asymptotic\" for a multivariate 'X' with ",
        "covariates.", call. = FALSE)
    }
    stopifnot(ncol(Z) < 2 | is.null(Z))

    # pool sized to the largest number of permutations any adaptive stage needs
    n_perm_pool <- ifelse(isTRUE(adaptive), sum(n_perm_adaptive), n_perm)
    X_star <- X_perm(X, Z, n_perm = n_perm_pool)

    if (adaptive == TRUE) {
      #### adaptive ----
      message(paste("Computing", n_perm_adaptive[1], "permutations..."))

      res0 <- pbapply::pblapply(seq_len(r),
        function(j) {
          cit_perm(
            Y = M[, j],
            X = X,
            Z = Z,
            X_star = X_star,
            n_perm = n_perm_adaptive[1],
            space_y = space_y,
            number_y = number_y)[c("score", "test_statistic")]
        },
        cl = par_clust)

      res  <- vapply(X = res0, FUN = function(u) {
        u$score
      },
      FUN.VALUE = numeric(1))
      stat <- vapply(X = res0, FUN = function(u) {
        u$test_statistic
      },
      FUN.VALUE = numeric(1))
      perm <- rep(n_perm_adaptive[1], r)
      used_perms <- n_perm_adaptive[1]

      k <- 2
      while (k <= length(n_perm_adaptive)) {
        index <- which(((res + 1) / (perm + 1)) < thresholds[k - 1])
        if (length(index) == 0L) { # no more genes
          break
        }

        message(paste("Computing", sum(n_perm_adaptive[k]), "additional permutations..."))

        if (parallel & .Platform$OS.type == "unix") {
          res_perm <- pbapply::pbsapply(seq_along(index),
            function(i) {
              cit_perm(Y = M[, index[i]],
                X = X,
                Z = Z,
                X_star = X_star[(used_perms + 1):(used_perms + n_perm_adaptive[k])],
                n_perm = n_perm_adaptive[k],
                space_y = space_y,
                number_y = number_y)$score
            },
            cl = par_clust, mc.preschedule = TRUE)
        } else {
          res_perm <- pbapply::pbsapply(seq_along(index),
            function(i) {
              cit_perm(Y = M[, index[i]],
                X = X,
                Z = Z,
                X_star = X_star[(used_perms + 1):(used_perms + n_perm_adaptive[k])],
                n_perm = n_perm_adaptive[k],
                space_y = space_y,
                number_y = number_y)$score
            },
            cl = par_clust)
        }

        res[index] <- res[index] + res_perm
        perm[index] <- perm[index] + rep(n_perm_adaptive[k], length(index))
        used_perms <- used_perms + n_perm_adaptive[k]
        k <- k + 1

      }

      pvals <- (res + 1) / (perm + 1)
      df <- data.frame(raw_pval = pvals,
        adj_pval = p.adjust(pvals, method = "BH"),
        test_statistic = stat)

      n_perm <- cumsum(n_perm_adaptive)

    } else {
      #### non-adaptive ----
      message(paste("Computing", n_perm, "permutations..."))

      if (parallel & .Platform$OS.type == "unix") {
        res <- do.call("rbind", pbapply::pblapply(seq_len(r),
          function(j) {
            cit_perm(Y = M[, j],
              X = X,
              Z = Z,
              X_star = X_star,
              n_perm = n_perm,
              space_y = space_y,
              number_y = number_y)
          },
          cl = par_clust, mc.preschedule = TRUE))
      } else {
        res <- do.call("rbind", pbapply::pblapply(seq_len(r),
          function(j) {
            cit_perm(Y = M[, j],
              X = X,
              Z = Z,
              X_star = X_star,
              n_perm = n_perm,
              space_y = space_y,
              number_y = number_y)
          },
          cl = par_clust))
      }

      # res <- as.vector(unlist(res))

      df <- data.frame(raw_pval = res$raw_pval,
        adj_pval = p.adjust(res$raw_pval, method = "BH"),
        test_statistic = res$test_statistic)

    }


  } else if (test == "asymptotic") {
    ## asymptotic ----
    n_perm <- NA
    # model.matrix / crossprod / solve / H depend only on (X, Z): build once,
    # not once per gene.
    design <- .cit_design(X, Z, n)
    res <- do.call("rbind", pbapply::pblapply(seq_len(r),
      function(j) {
        cit_asymp(M[, j], X, Z,
          space_y = space_y,
          number_y = number_y,
          design = design)
      },
      cl = par_clust)
    )
    df <- data.frame(raw_pval = res$raw_pval,
      adj_pval = p.adjust(res$raw_pval, method = "BH"),
      test_statistic = res$test_statistic)
  }

  if (parallel && .Platform$OS.type != "unix") {
    parallel::stopCluster(par_clust)
  }

  if (!is.null(M_colnames)) {
    rownames(df) <- M_colnames
  }

  output <- list(which_test = test,
    n_perm = n_perm,
    pvals = df)
  class(output) <- "cit_multi"
  return(output)


}

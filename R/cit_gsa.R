# Design quantities depending only on (X, Z): constant across genes and,
# when Z is absent, across permutations of X.
.cit_gsa_design <- function(X, Z = NULL, n) {
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

# One gene: observed statistic + the n_perm permuted statistics, indexed by
# permutation so that entry k across genes all come from the SAME X_star[[k]].
.gene_perm_stats <- function(Y, Z, X_star, n_perm, space_y, number_y, design) {

  n  <- length(Y)
  Y  <- as.numeric(Y)
  oY <- order(Y)
  y <- .cit_y_grid(Y, space_y, number_y)
  p  <- length(y)

  ij <- findInterval(y[-p], sort.int(Y))    # O(n log n); == sapply(sum(Y <= .))

  stat <- function(Hm) {
    b <- c(apply(Hm[, oY, drop = FALSE], MARGIN = 1, FUN = cumsum)[ij, , drop = FALSE]) / n
    return(sum(b^2) * n)
  }

  obs  <- stat(design$H)
  perm <- numeric(n_perm)

  for (k in seq_len(n_perm)) {
    mmk <- if (is.null(Z)) model.matrix(~., data = X_star[[k]])
    else            model.matrix(~., data = cbind(X_star[[k]], Z))
    Hk  <- if (is.null(Z)) n * tcrossprod(design$XtXinv_X, mmk)
    else            n * (solve(crossprod(mmk)) %*% t(mmk))[design$indexes_X, , drop = FALSE]
    perm[k] <- stat(Hk)
  }

  return(list(obs = obs, perm = perm))
}

# One gene set: sum the per-gene statistics, then count how often the
# shared-permutation null sum meets or exceeds the observed sum.
.gsa_perm_set <- function(M, genes, Z, X_star, n_perm, space_y, number_y, design,
                          set_index = NA) {

  measured <- intersect(colnames(M), genes)

  if (length(measured) < 1L) {
    warning("0 genes from geneset ", set_index, " observed in expression data")
    return(list(score = NA_integer_, obs = NA_real_))
  }

  if (length(measured) < length(genes)) {
    warning(" Some genes from geneset ", set_index, " are not observed in expression data")
  }

  pg  <- lapply(measured, function(g)
    .gene_perm_stats(M[, g], Z, X_star, n_perm, space_y, number_y, design))
  obs <- sum(vapply(pg, `[[`, numeric(1), "obs"))
  # column j = gene j's null over the shared pool, so row k = the gene-set
  # statistic under permutation k.
  perm_sum <- rowSums(vapply(pg, `[[`, numeric(n_perm), "perm"))

  return(list(score = sum(perm_sum >= obs), obs = obs))
}

#' Conditional independence test for gene set analysis
#'
#' @param M a \code{data.frame} or a \code{matrix} of size \code{n x r}
#' containing the different Y variables to test for conditional independence
#' with \code{X} adjusted on \code{Z}.
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
#' @param geneset a vector, a list, a gmt file format or a BiocSet object.
#' If the parameter is \itemize{
#'  \item a vector : corresponds to the gene name of the gene set, must be the same as those of the columns of the matrix \code{M}
#'  \item a list : each elements of the list are a gene set with the names of the genes, must be the same as those of the columns of the matrix \code{M}
#'  \item a gmt file format : the genes names of each genes set in the file, must be the same as those of the columns of the matrix \code{M}
#'  \item a BiocSet object : the genes names of each genes set in the object, must be the same as those of the columns of the matrix \code{M}
#' }
#'
#' @param test a character string indicating whether the \code{'asymptotic'} or
#' the \code{'permutation'} test is computed.
#' Default is \code{'asymptotic'}.
#'
#' @param n_perm the number of permutations. Default is \code{100}. Only used if
#' \code{test == 'permutation'}.
#'
#' @param adaptive a logical flag indicating whether adaptive permutations
#' should be performed. Default is \code{FALSE}. Only used if
#' \code{test == 'permutation'}.
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
#' \code{max(1L, paralleldetectCores(logical = FALSE) - 1L, na.rm = TRUE)}. If
#' \code{n_cpus = 1}, then sequential computations are used without any
#' parallelization.
#'
#' @param space_y a logical flag indicating whether the y thresholds are spaced out.
#' When \code{space_y} is \code{TRUE}, a regular sequence between the minimum and
#' the maximum of the observations is used. If \code{FALSE}, each unique
#' observed expression value is used as a distinct threshold. Default is \code{TRUE}.
#'
#' @param number_y an integer value indicating the number of y thresholds (and therefore
#' the number of regressions) to perform the test. Default is 10.
#'
#' @return A list with the following elements:\itemize{
#'   \item \code{which_test}: a character string carrying forward the value of
#'   the '\code{test}' argument indicating which test was performed (either
#'   'asymptotic' or 'permutation').
#'   \item \code{n_perm}: an integer carrying forward the value of the
#'   '\code{n_perm}' argument or '\code{n_perm_adaptive}' indicating the number of permutations performed
#'   (\code{NA} if asymptotic test was performed).
#'   \item \code{pvals}: computed p-values. A data frame with one row for
#'   each gene set, and with 2 columns: the first one '\code{raw_pval}' contains
#'   the raw p-values, the second one '\code{adj_pval}' contains the FDR adjusted p-values
#'   using Benjamini-Hochberg correction. When '\code{test == "asymptotic"}', a
#'   third column '\code{test_statistic}' contains the gene set test statistics.
#'   Gene sets with no gene observed in \code{M} yield a warning and \code{NA} in
#'   every column; gene sets only partially observed yield a warning and are
#'   tested on the measured genes alone.
#'   \item \code{type}: a character string equal to \code{"gsa"}, identifying
#'   the object as the result of a gene set analysis.
#' }
#'
#' @details The gene-set statistic is the sum of per-gene statistics. For the
#' permutation test, it is computed with each single permutation of X shared and
#' applied across all genes in a set (so inter-gene correlation is preserved).
#'
#' @seealso \code{\link{cit_multi}}, \code{\link{plot.cit_gsa}}
#'
#' @export
#'
#' @examples
#' set.seed(123)
#' n <- 500
#' r <- 200
#' Z1 <- rnorm(n)
#' Z2 <- rnorm(n) # rbinom(n, size=1, prob=0.5) + rnorm(n, sd=0.05)
#' X1 <- Z2 + rnorm(n, sd = 0.2)
#' X2 <- rnorm(n)
#' cor(X1, Z2)
#' Y <- replicate(r, Z2) + rnorm(n * r, 0, 0.5)
#' range(cor(Y, Z2))
#' range(cor(Y, X2))
#' res_asymp_unadj <- cit_gsa(M = data.frame(Y = Y),
#'   X = data.frame(X2 = X1),
#'   geneset = paste0("Y.", 1:50),
#'   test = "asymptotic", parallel = FALSE)
#' res_asymp_unadj$pvals
#'
#' # permutation test on a small example
#' set.seed(123)
#' n <- 60
#' X <- data.frame(X = as.factor(rbinom(n, size = 1, prob = 0.5)))
#' M <- matrix(rnorm(n * 8), nrow = n,
#'   dimnames = list(NULL, paste0("g", 1:8)))
#' geneset <- list(set1 = paste0("g", 1:4), set2 = paste0("g", 5:8))
#' res_perm <- cit_gsa(M = M, X = X, geneset = geneset,
#'   test = "permutation", n_perm = 50, parallel = FALSE)
#' res_perm$pvals
cit_gsa <- function(M,
                    X,
                    Z = NULL,
                    geneset,
                    test = c("asymptotic", "permutation"),
                    n_perm = 100,
                    n_perm_adaptive = c(n_perm, n_perm, n_perm * 3, n_perm * 5),
                    thresholds = c(0.1, 0.05, 0.01),
                    parallel = interactive(),
                    n_cpus = max(1L, detectCores(logical = FALSE) - 1L, na.rm = TRUE),
                    adaptive = FALSE,
                    space_y = TRUE,
                    number_y = 10) {
  # checks

  stopifnot(is.data.frame(M) | is.matrix(M))
  stopifnot(is.data.frame(X))
  stopifnot(is.data.frame(Z) | is.null(Z))
  stopifnot(is.logical(parallel))
  stopifnot(is.logical(adaptive))
  stopifnot(is.numeric(n_perm))
  stopifnot(inherits(geneset, "GSA.genesets") | inherits(geneset, "BiocSet") | is.character(geneset) | is.list(geneset))


  M_colnames <- colnames(M)


  if (anyNA(M)) {
    warning("'M' contains", sum(is.na(M)), "NA values. ",
      "\nCurrently they are ignored in the computations but ",
      "you should think carefully about where do those NA/NaN ",
      "come from...")
    M <- M[, colSums(is.na(M)) == 0]
  }

  r <- ncol(M)
  n <- nrow(M)
  stopifnot(nrow(X) == n)
  stopifnot(nrow(Z) == n | is.null(Z))

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
    if (test == "asymptotic") {
      n_cpus <- min(length(geneset), n_cpus)
    }
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

    if (!is.null(Z)) stopifnot(ncol(X) < 2)      # multi-column X supported only without Z
    stopifnot(ncol(Z) < 2 | is.null(Z))

    # normalise geneset
    if (inherits(geneset, "GSA.genesets")) {
      geneset <- geneset$genesets
    } else if (inherits(geneset, "BiocSet")) {
      if (!requireNamespace("BiocSet", quietly = TRUE)) {
        stop("Package 'BiocSet' is required for BiocSet input. Please install it from Bioconductor.")
      }
      geneset <- BiocSet::es_elementset(geneset)
      geneset <- lapply(unique(geneset$set), function(x) geneset[geneset$set == x, ]$element)
    } else if (is.vector(geneset) & !is.list(geneset)) {
      geneset <- list(geneset)
    }

    n_perm_pool <- if (isTRUE(adaptive)) sum(n_perm_adaptive) else n_perm
    X_star <- X_perm(X, Z, n_perm = n_perm_pool)

    design <- .cit_gsa_design(X, Z, n)

    if (isTRUE(adaptive)) {
      #### adaptive ----
      message(paste("Computing", n_perm_adaptive[1], "permutations..."))

      res0 <- pbapply::pblapply(seq_along(geneset), function(k)
        .gsa_perm_set(M, geneset[[k]], Z, X_star[seq_len(n_perm_adaptive[1])],
          n_perm_adaptive[1], space_y, number_y, design, set_index = k),
      cl = par_clust)

      score <- vapply(res0, `[[`, numeric(1), "score")
      obs   <- vapply(res0, `[[`, numeric(1), "obs")
      perm  <- rep(n_perm_adaptive[1], length(geneset))
      used_perms  <- n_perm_adaptive[1]

      k <- 2
      while (k <= length(n_perm_adaptive)) {
        index <- which(((score + 1) / (perm + 1)) < thresholds[k - 1])
        if (length(index) == 0) break
        message(paste("Computing", n_perm_adaptive[k], "additional permutations..."))

        slice <- X_star[(used_perms + 1):(used_perms + n_perm_adaptive[k])]   # disjoint
        res_k <- pbapply::pblapply(index, function(i)
          .gsa_perm_set(M, geneset[[i]], Z, slice,
            n_perm_adaptive[k], space_y, number_y, design, set_index = i),
        cl = par_clust)

        score[index] <- score[index] + vapply(res_k, `[[`, numeric(1), "score")
        perm[index]  <- perm[index] + n_perm_adaptive[k]
        used_perms <- used_perms + n_perm_adaptive[k]
        k <- k + 1
      }

      pvals <- (score + 1) / (perm + 1)
      df <- data.frame(raw_pval = pvals,
        adj_pval = p.adjust(pvals, method = "BH"),
        test_statistic = obs)
      n_perm <- cumsum(n_perm_adaptive)

    } else {
      #### non-adaptive ----
      message(paste("Computing", n_perm, "permutations..."))

      res <- pbapply::pblapply(seq_along(geneset), function(k)
        .gsa_perm_set(M, geneset[[k]], Z, X_star, n_perm,
          space_y, number_y, design, set_index = k),
      cl = par_clust)

      score <- vapply(res, `[[`, numeric(1), "score")
      pvals <- (score + 1) / (n_perm + 1)
      df <- data.frame(raw_pval = pvals,
        adj_pval = p.adjust(pvals, method = "BH"),
        test_statistic = vapply(res, `[[`, numeric(1), "obs"))
    }



  } else if (test == "asymptotic") {
    ## asymptotic ----
    n_perm <- NA



    # Data formatting in list format +  check column names
    if (inherits(geneset, "GSA.genesets")) {
      geneset <- geneset$genesets
    } else if (inherits(geneset, "BiocSet")) {
      if (!requireNamespace("BiocSet", quietly = TRUE)) {
        stop("Package 'BiocSet' is required for BiocSet input. Please install it from Bioconductor.")
      }
      geneset <- BiocSet::es_elementset(geneset)
      geneset <- lapply(X  = unique(geneset$set),
        FUN = function(x) {
          geneset[geneset$set == x, ]$element
        }
      )
    } else if (is.vector(geneset) & !is.list(geneset)) {
      geneset <- list(geneset)
    }

    # Initialisation for each gene set
    test_stat_list <- list()
    pval <- NA

    # design depends only on X and Z, and is not gene specific: compute it only once !
    n_Y_all <- nrow(M)
    H <- .cit_gsa_design(X, Z, n_Y_all)$H


    if (length(geneset) < 3) {
      op <- pbapply::pboptions(type = "none")
      on.exit(pbapply::pboptions(op), add = TRUE)
    }


    res <- pbapply::pblapply(seq_along(geneset), function(k) { # 1 -- each list of gene set ----

      # Initialisation for each gene in the gene set
      # test_stat_gs <- NULL
      prop_gs <- list()
      indi_pi_gs <- list()
      # ccdf_gs <- list()


      measured_genes <- intersect(M_colnames, geneset[[k]])

      if (length(measured_genes) < 1) { # check 1 : none genes of the current geneset are in M
        warning("0 genes from geneset ", k, " observed in expression data")
        pval <- NA
        test_stat_gs <- NA
        # test_stat_list <- NA
      } else {

        test_stat_gs <- numeric(length(measured_genes))

        if (length(measured_genes) < length(geneset[[k]])) { # check 2 : some genes of the current geneset are not in M
          warning(" Some genes from geneset ", k, " are not observed in expression data")
        }


        # code below if all the genes in the gene set k are in M, else pval + stat de test = NA

        for (i in seq_along(measured_genes)) { # 2 -- each genes in the gene set k ----

          Y <- M[, measured_genes[i]]
          oY <- order(Y)


          # 1) Test statistic computation ----
          y <- .cit_y_grid(Y, space_y, number_y)
          p <- length(y)

          index_jumps <- findInterval(y[-p], Y[oY])
          beta <- c(apply(X = H[, oY, drop = FALSE], MARGIN = 1, FUN = cumsum)[index_jumps, , drop = FALSE]) / n_Y_all # same number than thresholds
          test_stat <- sum(beta^2) * n_Y_all

          test_stat_gs[i] <- test_stat # test statistic for each genes in the gene set


          # 2) Pi computation ----
          indi_pi <- outer(Y, y[-p], "<=") * 1

          indi_pi_gs[[i]] <- indi_pi
          prop <- colMeans(indi_pi)
          prop_gs[[i]] <- prop # prop for each genes in the gene set


          # ccdf_gs[[i]] <- ccdf(Y=Y, X=X, Z=Z, method="OLS", fast=TRUE, space_y=space_y, number_y=number_y)


        }
        # ccdf_list[[k]] <- ccdf_gs
        # names(ccdf_list[[k]]) <- measured_genes # not geneset, if some genes are not in the data
        # utile comme le refait plus tard ???


        indi_pi_gs_tab <- do.call(cbind, indi_pi_gs)
        prop_gs_vec <- unlist(prop_gs)
        n_g_t <- length(prop_gs_vec)

        # 3) Sigma matrix creation ----
        n_gs_vec <- nrow(indi_pi_gs_tab)
        temp <- indi_pi_gs_tab - matrix(prop_gs_vec, nrow = n_gs_vec, ncol = n_g_t, byrow = TRUE)

        # `temp` is already centred, so crossprod(temp)/n_gs_vec IS the geneset covariance
        covmat <- crossprod(temp) / n_gs_vec

        ev_H   <- eigen(tcrossprod(H), symmetric = TRUE, only.values = TRUE)$values
        ev_cov <- eigen(covmat, symmetric = TRUE, only.values = TRUE)$values
        ev     <- as.vector(outer(ev_H, ev_cov)) / n

        pval <- survey::pchisqsum(sum(test_stat_gs), lower.tail = FALSE,
          df = rep(1, length(ev)),
          a = ev, method = "saddlepoint")
      }

      return(list("pval" = pval, "test_stat_gs" = test_stat_gs)) # ,"ccdf" = ccdf_list

    },
    cl = par_clust)

    pvals <- sapply(res, "[[", "pval")

    test_stat_list <- lapply(res, "[[", "test_stat_gs")

    df <- data.frame(raw_pval = pvals,
      adj_pval = p.adjust(pvals, method = "BH"),
      test_statistic = sapply(test_stat_list, sum))
  }

  if (parallel && .Platform$OS.type != "unix") {
    parallel::stopCluster(par_clust)
  }

  output <- list(which_test = test,
    n_perm = n_perm,
    pvals = df)

  class(output) <- "cit_gsa"
  output$type <- "gsa"
  return(output)

}

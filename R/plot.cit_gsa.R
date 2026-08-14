# =============================================================================
# Internal helpers for the gene-set CCDF plot
#
# These reuse the conventions established in plot_compare_ccdf.R: viridis for
# the levels of X, `goldenrod1` dashed for the reference curve, and
# `.ccdf_labs()` / `.ccdf_stack()` for the shared theme and A/B stacking.
# =============================================================================

# Evaluate a step function on an arbitrary grid.
#
# The CCDF of each gene is a right-continuous step function observed on that
# gene's own expression values. Genes in a set do not share a grid, so every
# curve is re-evaluated on a common one before they can be summarized. Grid
# Below a gene's first observation the value is 0 -- an empirical CDF is zero
# below its support, so the curve IS defined there. Returning NA instead would
# make the summary take a median over a changing subset of genes, which breaks
# monotonicity as each new gene enters the grid with the low start of its own
# CDF. Beyond the last observation the final value is carried forward -- the
# step-function extension the previous implementation did by hand.
.gsa_step_eval <- function(y_obs, v_obs, grid) {
  o   <- order(y_obs)
  yo  <- y_obs[o]
  vo  <- v_obs[o]
  idx <- findInterval(grid, yo)
  out <- rep(0, length(grid))
  ok  <- idx >= 1L
  out[ok] <- vo[idx[ok]]
  out
}

# Per-gene CCDF curves for one gene set.
#
# `by` selects the conditioning used for the returned curves:
#   "x"  -> CCDF given X only (panel A)
#   "xz" -> CCDF given X and Z, via a single interaction factor (panel B)
#   "z"  -> CCDF given Z only, i.e. the reference marginal on Z
# Passing a single factor to ccdf() makes the underlying OLS fit saturated, so
# each returned curve is an exact within-cell empirical CDF: monotonic and
# confined to [0, 1]. See ?plot_compare_ccdf for why that matters.
.gsa_curves <- function(M, genes, Xd, Zd, by, method, fast, space_y, number_y) {
  grp <- switch(by,
    "x"  = Xd,
    "z"  = Zd,
    "xz" = interaction(Xd, Zd, sep = "|", drop = TRUE))

  do.call(rbind, lapply(genes, function(g) {
    res <- ccdf(as.numeric(M[, g]), data.frame(X = grp), NULL,
      method, fast, space_y, number_y)
    lab <- as.character(res$x)
    if (by == "xz") {
      sp <- do.call(rbind, strsplit(lab, "|", fixed = TRUE))
      data.frame(gene = g, y = res$y, ccdf = res$ccdf, cdf = res$cdf,
        x = sp[, 1], z = sp[, 2], stringsAsFactors = FALSE)
    } else {
      data.frame(gene = g, y = res$y, ccdf = res$ccdf, cdf = res$cdf,
        x = lab, z = NA_character_, stringsAsFactors = FALSE)
    }
  }))
}

# Pointwise summary across genes, on a shared grid.
#
# Each gene's step function is evaluated on `grid`, then summarized within every
# (x, z) cell. Because every input curve is monotonic (see .gsa_curves), a
# pointwise median of them is monotonic too -- so the summary curve needs no
# monotonicity correction, unlike the previous implementation.
.gsa_summary <- function(df, grid, value = "ccdf", fun = stats::median) {
  df$z2  <- ifelse(is.na(df$z), "__none__", df$z)
  cells  <- unique(df[, c("x", "z2")])
  do.call(rbind, lapply(seq_len(nrow(cells)), function(i) {
    sub <- df[df$x == cells$x[i] & df$z2 == cells$z2[i], , drop = FALSE]
    mat <- vapply(split(sub, sub$gene),
      function(d) .gsa_step_eval(d$y, d[[value]], grid),
      numeric(length(grid)))
    if (is.null(dim(mat))) mat <- matrix(mat, nrow = length(grid))
    data.frame(y = grid,
      value = apply(mat, 1, function(v) fun(v[!is.na(v)])),
      x = cells$x[i],
      z = if (cells$z2[i] == "__none__") NA_character_ else cells$z2[i],
      stringsAsFactors = FALSE)
  }))
}

# Build one gene-set panel: transparent per-gene curves + bold summaries.
.gsa_panel <- function(curves, summ, ref, lv, keys, y_lab, alpha, linewidth,
                       legend_name = "CCDF", facet = FALSE,
                       facet_labeller = NULL) {
  vx      <- viridis(n = length(lv) + 1)[seq_len(length(lv))]
  ref_lbl <- if (facet) "Marginal on X" else "CDF"
  curves$x <- factor(curves$x, levels = lv)
  summ$x   <- factor(summ$x,   levels = lv)

  p <- ggplot(mapping = aes(x = .data$y)) +
    # individual genes: semi-transparent, one line per gene per level of X
    geom_step(data = curves,
      aes(y = .data$ccdf, color = .data$x,
        group = interaction(.data$gene, .data$x)),
      linewidth = 0.25, alpha = alpha, show.legend = FALSE) +
    # gene-set summary: same color, opaque and thicker
    geom_step(data = summ,
      aes(y = .data$value, color = .data$x), linewidth = linewidth) +
    # reference curve, matching plot_compare_ccdf's gold dashed convention
    geom_step(data = ref,
      aes(y = .data$value, color = ref_lbl),
      linewidth = linewidth, linetype = "dashed") +
    scale_color_manual(
      name   = legend_name,
      values = stats::setNames(c(vx, "goldenrod1"), c(lv, ref_lbl)),
      breaks = c(lv, ref_lbl),
      labels = c(keys, ref_lbl),
      guide  = guide_legend(override.aes = list(
        linetype = c(rep("solid", length(lv)), "dashed"),
        shape    = rep(NA, length(lv) + 1)))) +
    .ccdf_labs(y_lab)  +
    theme_classic()

  if (facet) {
    p <- p + facet_wrap(~z, ncol = 2, scales = "fixed", labeller = facet_labeller)
  }
  return(p)
}


# =============================================================================
# Main S3 method
# =============================================================================

#' Plot the conditional CDFs of every gene in a gene set
#'
#' Draws, for one gene set, the conditional CDF of each gene given \code{X}
#' (and \code{Z}) as a semi-transparent step function, overlaid with a bold
#' gene-set summary curve per level of \code{X}. Layout, colors, reference
#' curve and faceting follow \code{\link{plot_compare_ccdf}}, of which this is
#' the many-genes counterpart.
#'
#' @param x an object of class \code{cit_gsa}, as returned by
#' \code{\link{cit_gsa}}.
#'
#' @param M a numeric matrix or data frame of size \code{n x r} containing the
#' preprocessed expressions, with gene names as column names. Required:
#' \code{cit_gsa()} does not retain the data it was called on.
#'
#' @param X a data frame whose first column is the variable of interest.
#'
#' @param Z a data frame whose first column is the covariate, or \code{NULL}.
#'
#' @param geneset the gene set to draw: a character vector of gene names, or a
#' named list of such vectors, in which case \code{which_set} picks one.
#'
#' @param which_set index or name selecting the gene set when \code{geneset} is
#' a list. Default is the first.
#'
#' @param discretize a logical flag, following the same logic as in
#' \code{\link{plot_compare_ccdf}}: continuous variables are quartile-binned and
#' combined into a single interaction factor, so every curve drawn is an exact
#' empirical CDF. Because a pointwise median of monotonic curves is itself
#' monotonic, the summary curve then needs no monotonicity correction. Default
#' is \code{FALSE} when \code{X} (and \code{Z}) are already factors, and
#' \code{TRUE} otherwise.
#'
#' @param summary_fun the function used to summarize across genes at each
#' expression value. Default is \code{\link[stats]{median}}.
#'
#' @param alpha opacity of the individual gene curves. Default \code{0.25}.
#'
#' @param linewidth width of the summary and reference curves. Default
#' \code{0.9}.
#'
#' @param n_grid number of points on the shared grid used to summarize across
#' genes. Genes in a set do not share expression values, so all curves are
#' re-evaluated on a common grid first. Default \code{200}.
#'
#' @param method,fast,space_y,number_y passed to \code{\link{ccdf}}.
#'
#' @param probs,bin_labels passed to the quartile binning, exactly as in
#' \code{\link{plot_compare_ccdf}}.
#'
#' @param ... further arguments to be passed.
#'
#' @return a \code{\link[ggplot2]{ggplot}} object; when \code{Z} is supplied, a
#' \code{\link[patchwork]{patchwork}} composition stacking the CCDF given
#' \code{X} alone (panel A) above the CCDF given \code{X} and \code{Z}
#' (panel B, faceted by \code{Z}).
#'
#' @importFrom viridisLite viridis
#' @importFrom stats median
#' @import ggplot2 patchwork
#'
#' @seealso \code{\link{plot_compare_ccdf}} for the single-gene version.
#'
#' @export
#'
#' @examples
#' set.seed(123)
#' n <- 60
#' X <- data.frame(X = as.factor(rbinom(n, size = 1, prob = 0.5)))
#'
#' # 20 genes: two sets of 10. Only the genes of set1 depend on X.
#' M <- matrix(rnorm(n * 20), nrow = n,
#'   dimnames = list(NULL, paste0("g", 1:20)))
#' M[, 1:10] <- M[, 1:10] + 1.5 * (as.numeric(X$X) - 1)
#' geneset <- list(set1 = paste0("g", 1:10), set2 = paste0("g", 11:20))
#'
#' res <- cit_gsa(M = M, X = X, geneset = geneset,
#'   test = "asymptotic", parallel = FALSE)
#'
#' # set1: the two summary curves separate, and the individual genes with them
#' plot(res, M = M, X = X, geneset = geneset, which_set = "set1")
#'
#' \donttest{
#' # set2 is null: the two summary curves stay close to each other and to the
#' # reference CDF
#' plot(res, M = M, X = X, geneset = geneset, which_set = "set2")
#'
#' # with a continuous covariate: quartile-binned and faceted by Z
#' Z <- data.frame(Z = rnorm(n))
#' plot(res, M = M, X = X, Z = Z, geneset = geneset, which_set = 1)
#' plot(res, M = M, X = X, Z = Z, geneset = geneset, which_set = 2)
#' }
#'
plot.cit_gsa <- function(x, M, X, Z = NULL, geneset, which_set = 1,
                         method = c("OLS", "logistic"), fast = TRUE,
                         space_y = FALSE, number_y = 20,
                         discretize = !is.factor(X[, 1]) ||
                           (!is.null(Z) && !is.factor(Z[, 1])),
                         probs = c(0, .25, .5, .75, 1),
                         bin_labels = c("Q1", "Q2", "Q3", "Q4"),
                         summary_fun = stats::median, alpha = 0.25,
                         linewidth = 0.9, n_grid = 200, ...) {

  stopifnot(inherits(x, "cit_gsa"))
  if (missing(M) || missing(X) || missing(geneset))
    stop("'M', 'X' and 'geneset' are required: a cit_gsa object does not ",
      "retain the data it was computed from.")

  ## ---- resolve the gene set ----------------------------------------------
  gs <- if (is.list(geneset)) {
    geneset[[which_set]]
  } else {
    geneset
  }
  # name of the selected set, used as the plot title when there is one
  set_name <- if (is.list(geneset)) {
    nms <- names(geneset)
    if (is.character(which_set)) {
      which_set
    } else if (!is.null(nms)) {
      nms[which_set]
    } else {
      NULL
    }
  } else {
    NULL
  }

  genes <- intersect(as.character(gs), colnames(M))
  if (length(genes) < 1) {
    stop("none of the genes in the selected set are present in the columns of 'M'.")
  }

  y_lab <- "Gene expression"
  Xd <- if (isTRUE(discretize)) .ccdf_qbin(X[, 1], probs, bin_labels) else X[, 1]
  Xd <- as.factor(Xd)
  lv <- levels(Xd)
  x_keys <- paste0("CCDF | ", .ccdf_bin_labels(lv, X[, 1], "X"))


  grid_from <- function(df) seq(min(df$y, na.rm = TRUE),
    max(df$y, na.rm = TRUE), length.out = n_grid)

  ## ---- panel A: given X, ignoring Z --------------------------------------
  cur_x <- .gsa_curves(M, genes, Xd, NULL, "x", method, fast, space_y, number_y)
  gridA <- grid_from(cur_x)
  sumA  <- .gsa_summary(cur_x, gridA, "ccdf", summary_fun)
  ref_x <- cur_x
  ref_x$x <- "all"
  refA  <- .gsa_summary(ref_x, gridA, "cdf", summary_fun)
  p_top <- .gsa_panel(cur_x, sumA, refA, lv, x_keys, y_lab, alpha, linewidth,
    legend_name = NULL, facet = FALSE)

  if (is.null(Z)) {
    if (!is.null(set_name)) {
      p_top <- p_top +
        ggtitle(set_name) +
        theme(plot.title = element_text(hjust = 0.5))
    }
    return(p_top)
  }
  p_top <- p_top + ggtitle("Marginal on Z")


  ## ---- panel B: given X and Z, faceted by Z ------------------------------
  Zd <- if (isTRUE(discretize)) .ccdf_qbin(Z[, 1], probs, bin_labels) else Z[, 1]
  Zd <- as.factor(Zd)
  z_lab <- colnames(Z)
  z_facet_labeller <- as_labeller(stats::setNames(paste("|",
    .ccdf_bin_labels(levels(Zd), Z[, 1], z_lab)), levels(Zd)))

  cur_xz <- .gsa_curves(M, genes, Xd, Zd, "xz", method, fast, space_y, number_y)
  gridB  <- grid_from(cur_xz)
  sumB   <- .gsa_summary(cur_xz, gridB, "ccdf", summary_fun)

  ## reference: CCDF marginal on Z, one curve per facet
  cur_z <- .gsa_curves(M, genes, Xd, Zd, "z", method, fast, space_y, number_y)
  refB  <- .gsa_summary(cur_z, gridB, "ccdf", summary_fun)
  refB$z <- refB$x
  refB$x <- NA_character_

  p_bottom <- .gsa_panel(cur_xz, sumB, refB, lv, x_keys, y_lab, alpha, linewidth,
    legend_name = "CCDF", facet = TRUE,
    facet_labeller = z_facet_labeller)

  out <- .ccdf_stack(p_top, p_bottom)
  if (!is.null(set_name)) {
    out <- out +
      plot_annotation(title = set_name,
        theme = theme(plot.title = element_text(hjust = 0.5)))
  }

  return(out)
}

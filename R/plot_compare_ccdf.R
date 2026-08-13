# =============================================================================
# Internal helpers -- not exported (leading "." keeps them out of NAMESPACE)
# =============================================================================

# ---------------------------------------------------------------- helpers ----
.ccdf_labs <- function(y_lab) {
  list(xlab(y_lab), ylab("Cumulative probability"), theme_bw(),
    theme(plot.title = element_text(hjust = 0.5)))
}

# One "series" = one drawn quantity. `labels` are the colour-group values,
# `keys` the text shown in the legend. A mapped series expands over levels of x.
.ccdf_serie <- function(ycol, labels, colours, geom, keys = labels,
                        mapped = FALSE, linetype = "solid", shape = 16,
                        linewidth = 0.5, size = 0.5) {
  list(ycol = ycol, labels = labels, keys = keys, colours = colours,
    geom = geom, mapped = mapped,
    linetype = if (geom == "point") "blank" else linetype,
    shape    = if (geom == "point") shape    else NA,
    linewidth = linewidth, size = size)
}

# Every reference curve (the CDF, or the CCDF marginal on Z) is drawn in this
# colour, dotted when it is a line, so it reads the same way in every panel.

.ccdf_ref <- function(ycol, label, geom = "step", ...) {
  .ccdf_serie(ycol, label, "goldenrod1", geom, linetype = "dashed", ...)
}

.ccdf_layer <- function(s) {
  if (s$mapped) {
    if (s$geom == "step")
      geom_step(aes(y = .data[[s$ycol]], color = .data$x),
        linewidth = s$linewidth, linetype = s$linetype)
    else
      geom_point(aes(y = .data[[s$ycol]], color = .data$x),
        shape = s$shape, size = s$size)
  } else {
    lbl <- s$labels[1]
    if (s$geom == "step")
      geom_step(aes(y = .data[[s$ycol]], color = lbl),
        linewidth = s$linewidth, linetype = s$linetype)
    else
      geom_point(aes(y = .data[[s$ycol]], color = lbl),
        shape = s$shape, size = s$size)
  }
}

# The scale is derived from the SAME spec as the layers, so legend colours,
# linetypes and shapes cannot drift away from what is drawn.
.ccdf_scale <- function(series, name = "") {
  pull <- function(f) unlist(lapply(series, f), use.names = FALSE)
  labs <- pull(function(s) s$labels)
  scale_color_manual(
    name   = name,
    values = stats::setNames(pull(function(s) s$colours), labs),
    breaks = labs,
    labels = pull(function(s) s$keys),
    guide  = guide_legend(override.aes = list(
      linetype = pull(function(s) rep(s$linetype, length(s$labels))),
      shape    = pull(function(s) rep(s$shape,    length(s$labels))))))
}

.ccdf_panel <- function(df, series, y_lab, title = NULL, legend_name = "") {
  p <- ggplot(df, aes(x = .data$y)) + lapply(series, .ccdf_layer) +
    .ccdf_scale(series, legend_name) + .ccdf_labs(y_lab)
  if (!is.null(title)) p <- p + ggtitle(title)
  p
}

.ccdf_stack <- function(top, bottom) {
  # The spacer reproduces the 75%-width top panel of the previous cowplot
  # layout; patchwork does not tag spacers, so the panels are still A and B.
  ((top | plot_spacer()) + plot_layout(widths = c(3, 1))) /
    bottom +
    plot_layout(heights = c(1, 3)) +
    plot_annotation(tag_levels = "A",
      theme = theme(plot.tag = element_text(size = 15)))
}

# Quartile-bin a continuous variable for the `discretize` display mode; a
# factor is returned unchanged. Ties can collapse quartile breaks, so
# unique() is applied to the break points -- this yields fewer than the
# requested bins rather than erroring.
.ccdf_qbin <- function(v, probs = c(0, .25, .5, .75, 1),
                       bin_labels = c("Q1", "Q2", "Q3", "Q4")) {
  if (is.factor(v)) return(v)
  brks <- unique(stats::quantile(v, probs, na.rm = TRUE))
  cut(v, breaks = brks, include.lowest = TRUE,
    labels = bin_labels[seq_len(length(brks) - 1)])
}

# Legend/facet text for a (possibly quartile-binned) variable.
#
# `lv` are the levels actually used for grouping (eg "Q1".."Q4" when binned),
# so they are replaced with the interval each bin actually spans, reusing the
# bin names themselves as the interior breakpoints: "X<Q1", "Q1\u2264X<Q2",
# ..., "Q(n-1)\u2264X".
.ccdf_bin_labels <- function(lv, orig, var_name) {
  if (is.factor(orig) || length(lv) < 2) return(paste0(var_name, "=", lv))
  cuts <- lv[-length(lv)]
  mid  <- if (length(cuts) > 1)
    paste0(cuts[-length(cuts)], "\u2264", var_name, "<", cuts[-1]) else character(0)
  return(c(paste0(var_name, "<", cuts[1]), mid, paste0(cuts[length(cuts)], "\u2264", var_name)))
}

# Panels for `discretize = TRUE`: X (and Z, if present) are quartile-binned
# (factors pass through unchanged), then combined into a SINGLE factor via
# interaction() before calling ccdf(). A saturated single-factor OLS fit
# returns exact within-cell proportions, so every curve is a genuine
# empirical CDF: monotonic and confined to [0, 1] by construction. This also
# repairs the rare case-C artifact (both X and Z already factors) -- see the
# roxygen note on plot_compare_ccdf().
.ccdf_discretize_panels <- function(Y, X, Z, method, fast, space_y, number_y,
                                    probs, bin_labels) {
  yv    <- as.numeric(Y[, 1])
  y_lab <- colnames(Y)
  Xd    <- .ccdf_qbin(X[, 1], probs, bin_labels)
  lv    <- levels(Xd)
  l_X <- length(lv)
  vx    <- viridis(n = l_X + 1)[seq_len(l_X)]

  if (is.null(Z)) {
    res <- ccdf(yv, data.frame(X = Xd), NULL, method, fast, space_y, number_y)
    df  <- data.frame(y = res$y, x = res$x, cdf = res$cdf, ccdf = res$ccdf)
    series <- list(.ccdf_ref("cdf", "CDF"),
      .ccdf_serie("ccdf", lv, vx, "step",
        keys = paste0("CCDF | ", .ccdf_bin_labels(lv, X[, 1], "X")), mapped = TRUE))
    return(list(.ccdf_panel(df, series, y_lab)))
  }

  Zd    <- .ccdf_qbin(Z[, 1], probs, bin_labels)
  z_lab <- colnames(Z)
  XZ    <- interaction(Xd, Zd, sep = "|", drop = TRUE)
  res   <- ccdf(yv, data.frame(X = XZ), NULL, method, fast, space_y, number_y)
  sp    <- do.call(rbind, strsplit(as.character(res$x), "|", fixed = TRUE))
  df    <- data.frame(y = res$y, ccdf = res$ccdf,
    x = factor(sp[, 1], levels = lv),
    z = factor(sp[, 2], levels = levels(Zd)))

  res_X <- ccdf(yv, data.frame(X = Xd), NULL, method, fast, space_y, number_y)
  df_X  <- data.frame(y = res_X$y, x = res_X$x, cdf = res_X$cdf, ccdf = res_X$ccdf)
  top_series <- list(.ccdf_serie("ccdf", lv, vx, "step",
    keys = paste0("CCDF | ", .ccdf_bin_labels(lv, X[, 1], "X")), mapped = TRUE),
  .ccdf_ref("cdf", "CDF"))
  p_top <- .ccdf_panel(df_X, top_series, y_lab) +
    ggtitle("Marginal on Z") +
    theme(plot.title = element_text(hjust = 0))

  res_Z <- ccdf(yv, data.frame(X = Zd), NULL, method, fast, space_y, number_y)
  df_Z  <- data.frame(y = res_Z$y, ccdf = res_Z$ccdf,
    z = factor(res_Z$x, levels = levels(Zd)))

  ref_lbl <- "Marginal on X"
  p_bottom <- ggplot(df, aes(x = .data$y)) +
    geom_step(aes(y = .data$ccdf, color = .data$x), linewidth = 0.5) +
    geom_step(data = df_Z, aes(y = .data$ccdf, color = ref_lbl),
      linewidth = 0.5, linetype = "dotted") +
    scale_color_manual(
      name   = "CCDF",
      values = stats::setNames(c(vx, "goldenrod1"), c(lv, ref_lbl)),
      breaks = c(lv, ref_lbl),
      labels = c(paste0("CCDF | ", .ccdf_bin_labels(lv, X[, 1], "X")), ref_lbl),
      guide  = guide_legend(override.aes = list(
        linetype = c(rep("solid", l_X), "dotted"),
        shape    = rep(NA, l_X + 1)))) +
    .ccdf_labs(y_lab) +
    facet_wrap(~z, ncol = 2, scales = "fixed",
      labeller = as_labeller(stats::setNames(
        paste("|", .ccdf_bin_labels(levels(Zd), Z[, 1], z_lab)),
        levels(Zd))))
  list(p_top, p_bottom)
}

# --------------------------------------------------------------- panels ------
.ccdf_panels <- function(Y, X, Z = NULL, method = c("OLS", "logistic"),
                         fast = TRUE, space_y = FALSE,
                         number_y = length(unique(Y[, 1])),
                         discretize = FALSE,
                         probs = c(0, .25, .5, .75, 1),
                         bin_labels = c("Q1", "Q2", "Q3", "Q4")) {
  if (isTRUE(discretize))
    return(.ccdf_discretize_panels(Y, X, Z, method, fast, space_y, number_y,
      probs, bin_labels))

  res   <- ccdf(as.numeric(Y[, 1]), X, Z, method, fast, space_y, number_y)
  y_lab <- colnames(Y)
  x_fac <- is.factor(X[, 1])
  lv    <- if (x_fac) levels(X[, 1]) else character(0)   # NOT unique(): see notes
  l_X   <- length(lv)
  vx    <- if (x_fac) viridis(n = l_X + 1) else viridis(n = 3)

  if (is.null(Z)) {
    df <- data.frame(y = res$y, x = res$x, cdf = res$cdf, ccdf = res$ccdf)
    series <- if (x_fac)
      list(.ccdf_ref("cdf", "CDF"),
        .ccdf_serie("ccdf", lv, vx[seq_len(l_X)], "step",
          keys = paste0("CCDF | X=", lv), mapped = TRUE))
    else
      list(.ccdf_ref("cdf", "CDF"),
        .ccdf_serie("ccdf", "CCDF", vx[1], "point"))
    return(list(.ccdf_panel(df, series, y_lab)))
  }

  df    <- data.frame(y = res$y, x = res$x, z = res$z, cdf = res$cdf,
    ccdf_nox = res$ccdf_nox, ccdf_x = res$ccdf_x)
  res_X <- ccdf(as.numeric(Y[, 1]), X, Z = NULL, method, fast, space_y, number_y)
  df_X  <- data.frame(y = res_X$y, x = res_X$x, cdf = res_X$cdf, ccdf = res_X$ccdf)
  z_fac <- is.factor(Z[, 1])
  z_lab <- colnames(Z)

  ## Z and X both continuous: everything fits on a single panel
  if (!z_fac && !x_fac) {
    series <- list(
      .ccdf_ref("cdf",      "CDF", "step", linewidth = 0.7),
      .ccdf_serie("ccdf_x",   "CCDF | Z, X",   vx[1],  "point"),
      .ccdf_serie("ccdf_nox", "CCDF | Z", vx[2], geom = "point", shape = 2))
    return(list(.ccdf_panel(df, series, y_lab)))
  }

  ## panel A -- marginal on Z
  if (x_fac) {
    top_series <- list(.ccdf_serie("ccdf", lv, vx[seq_len(l_X)], "step",
      keys = paste0("CCDF | X=", lv), mapped = TRUE),
    .ccdf_ref("cdf", "CDF"))
  } else {
    top_series <- list(.ccdf_serie("ccdf", "CCDF", vx[1], "point"),
      .ccdf_ref("cdf", "CDF"))
  }
  p_top <- .ccdf_panel(df_X, top_series, y_lab) +
    ggtitle("Marginal on Z") +
    theme(plot.title = element_text(hjust = 0))

  ## panel B -- given X and Z, vs marginal on Z
  bottom_series <- if (x_fac)
    list(.ccdf_serie("ccdf_x", lv, vx[seq_len(l_X)],
      if (z_fac) "step" else "point",
      keys = paste0("CCDF | X=", lv), mapped = TRUE),
    .ccdf_ref("ccdf_nox", "Marginal on X",
      geom = if (z_fac) "step" else "point", shape = 2))
  else
    list(.ccdf_serie("ccdf_x", "Given X and Z", vx[1], "point"),
      .ccdf_ref("ccdf_nox", "Marginal on X", geom = "point", shape = 2))

  p_bottom <- .ccdf_panel(df, bottom_series, y_lab,
    legend_name = "CCDF")
  if (z_fac) p_bottom <- p_bottom +
    facet_wrap(~z, ncol = if (x_fac) 3 else NULL,
      nrow = if (x_fac) NULL else 3, scales = "fixed",
      labeller = as_labeller(function(v) paste0("| ", z_lab, "=", v)))
  list(p_top, p_bottom)
}

# =============================================================================
# Main function
# =============================================================================

#' Function to plot the CCDF according to the type of X and Z
#'
#' @param Y a data frame whose first column contains the preprocessed
#' expressions from \code{n} samples (or cells). Its column name is used as the
#' x-axis label.
#'
#' @param X a data frame whose first column is a numeric or factor vector of
#' size \code{n} containing the variable to be tested (the condition to be
#' tested).
#'
#' @param Z a data frame whose first column is a numeric or factor vector of
#' size \code{n} containing the covariate. Multiple variables are not allowed.
#'
#' @param method a character string indicating which method to use to
#' compute the CCDF, either \code{'OLS'} or \code{'logistic'}.
#' Default is \code{'OLS'} for computational speed.
#'
#' @param fast a logical flag indicating whether the fast implementation of
#' logistic regression should be used. Only if \code{method == 'logistic'}.
#' Default is \code{TRUE}.
#'
#' @param space_y a logical flag indicating whether the y thresholds are spaced.
#' When \code{space_y} is \code{TRUE}, a regular sequence between the minimum and
#' the maximum of the observations is used. Default is \code{FALSE}.
#'
#' @param number_y an integer value indicating the number of y thresholds (and
#' therefore the number of regressions) used to compute the CCDF. Default is
#' \code{length(unique(Y[, 1]))}, i.e. one threshold per distinct observed value.
#'
#' @param discretize a logical flag.
#' When \code{TRUE}, any continuous variable
#' among \code{X} and \code{Z} is cut at \code{probs} into ordered bins. If
#' \code{Z} is not \code{NULL}, \code{X} and \code{Z} are combined into a single
#' interaction factor before calling \code{\link{ccdf}}. Default is \code{FALSE}
#' when \code{X} and \code{Z} are already factors, and \code{TRUE} otherwise.
#'
#' @section Discretization of continuous \code{X} and \code{Z}:
#'
#' \code{ccdf()} fits a separate regression at every \code{y} threshold. With no
#' constraint linking the the fits across thresholds, fitted
#' \eqn{P(Y <= y | X, Z)} are not necessarily monotonic whenever \code{X} and/or
#' \code{Z} are continuous variables: at each threshold the observed values for
#' \code{X} and \code{Z} varies, making the empirical conditioning different.
#' Such \code{ccdf} computations across varying covariate values at varying
#' thresholds do not carry any monotonicity guarantee, and become hard to interpret graphically.
#' For this reason, we provide the option to discretize \code{X} and
#' \code{Z} for graphical representation, in order to ease the interpretation.
#'
#' @section A note on when \code{X} and \code{Z} are both factors:
#'
#' When \code{X} and \code{Z} are already both factors, \code{ccdf()} fits
#' an *additive* model without an interaction term. With more than a handful of
#' levels in either variable this means that the fitted \code{ccdf_x} can sometimes
#' display some mild non-monotonicity (in practice this happens only with enough
#' levels on both sides to matter). Set \code{discretize = TRUE} to force the
#' saturated interaction encoding.
#'
#' @param probs breakpoints (as quantile probabilities) used to bin a continuous
#' \code{X} or \code{Z} when \code{discretize = TRUE}. Default is quartiles.
#' Ignored for variables that are already factors.
#'
#' @param bin_labels labels for the bins produced by \code{probs}. Default is
#' \code{c("Q1", "Q2", "Q3", "Q4")}. Ignored for variables that are already
#' factors, and truncated if \code{probs} produces fewer bins (ties can collapse
#' quantile breakpoints).
#'
#' @import ggplot2 patchwork stats
#' @importFrom viridisLite viridis
#'
#' @return a \code{\link[ggplot2]{ggplot}} object. When \code{Z} is supplied and
#' at least one of \code{X} and \code{Z} is a factor, the returned object is a
#' \code{\link[patchwork]{patchwork}} composition stacking the CCDF marginal on
#' \code{Z} (panel A) above the CCDF given \code{X} and \code{Z} (panel B).
#'
#' @export
#'
#' @examples
#' set.seed(123)
#' n <- 40
#' Y  <- data.frame(Y = rnorm(n))
#' Xf <- data.frame(X = as.factor(rbinom(n, size = 1, prob = 0.5)))
#' Xc <- data.frame(X = rnorm(n))
#' Zf <- data.frame(Z = as.factor(rbinom(n, size = 1, prob = 0.5)))
#' Zc <- data.frame(Z = rnorm(n))
#'
#' # 1. Z absent, X factor         -- CDF plus one CCDF step per level of X
#' plot_compare_ccdf(Y, Xf)
#'
#' # 2. Z absent, X continuous     -- CDF step plus CCDF points
#' plot_compare_ccdf(Y, Xc)
#'
#' # 3. Z factor, X factor         -- panel B faceted by Z, steps
#' plot_compare_ccdf(Y, Xf, Zf)
#'
#' # 4. Z factor, X continuous     -- panel B faceted by Z, points
#' plot_compare_ccdf(Y, Xc, Zf)
#'
#' # 5. Z continuous, X factor     -- panel B not faceted
#' plot_compare_ccdf(Y, Xf, Zc)
#'
#' # 6. Z continuous, X continuous -- a single panel, CDF plus both CCDFs
#' plot_compare_ccdf(Y, Xc, Zc)
#'
#' \donttest{
#' # a factor with more than two levels gets one colour per level
#' X3 <- data.frame(X = as.factor(sample(0:2, n, replace = TRUE)))
#' plot_compare_ccdf(Y, X3, Zf)
#'
#' # the logistic method works for every case above
#' plot_compare_ccdf(Y, Xf, Zf, method = "logistic")
#'
#' # force the interaction encoding even when X and Z are already factors
#' # (removes the rare additive-model artifact noted above)
#' plot_compare_ccdf(Y, Xf, Zf, discretize = TRUE)
#' }
#'
plot_compare_ccdf <- function(Y, X, Z = NULL, method = c("OLS", "logistic"),
                              fast = TRUE, space_y = FALSE,
                              number_y = length(unique(Y[, 1])),
                              discretize = !is.factor(X[, 1]) ||
                                (!is.null(Z) && !is.factor(Z[, 1])),
                              probs = c(0, .25, .5, .75, 1),
                              bin_labels = c("Q1", "Q2", "Q3", "Q4")) {
  panels <- .ccdf_panels(Y, X, Z, method, fast, space_y, number_y,
    discretize, probs, bin_labels)
  if (length(panels) == 1L) {
    return(panels[[1]])
  } else {
    return(.ccdf_stack(panels[[1]], panels[[2]]))
  }
}

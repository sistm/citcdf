#' citcdf: Conditional Independence Testing with Cumulative Distribution Functions
#'
#' Distribution-free conditional independence testing built on estimates of the
#' conditional cumulative distribution function (CCDF).
#'
#' @section Main functions:
#' \itemize{
#'   \item \code{\link{cit_multi}}: gene-wise testing across many outcomes.
#'   \item \code{\link{cit_gsa}}: gene-set analysis.
#'   \item \code{\link{cit_asymp}} / \code{\link{cit_perm}}: the single-outcome
#'         asymptotic and permutation tests.
#'   \item \code{\link{ccdf}}: the CCDF estimator the tests are built on.
#'   \item \code{\link{plot_compare_ccdf}}: diagnostic CCDF plots.
#' }
#'
#' @section Note on defaults:
#' \code{cit_multi()} and \code{cit_gsa()} use \code{space_y = TRUE} with
#' \code{number_y = 10} (for computational speed),
#' while \code{ccdf()}, \code{cit_asymp()} and \code{cit_perm()} default
#' to \code{space_y = FALSE}, i.e. every distinct observed value is a threshold.
#'
#' @references Gauthier M, Agniel D, Thiébaut R & Hejblum BP (2021).
#' Distribution-free complex hypothesis testing for single-cell RNA-seq
#' differential expression analysis, \emph{bioRxiv} 445165.
#' \doi{10.1101/2021.05.21.445165}

#' @keywords internal
"_PACKAGE"

## usethis namespace: start
## usethis namespace: end
NULL

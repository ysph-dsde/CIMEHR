#' Summary-Statistic Regression for Longitudinal Data
#'
#' @description
#' Reduces longitudinal data to a single summary statistic per subject
#' (e.g., mean, median, min, max of the outcome), then fits a
#' cross-sectional linear model on the summarised data. This provides a
#' simple baseline analysis that ignores the within-subject correlation
#' structure and visit process entirely.
#'
#' @param data A \code{data.frame} in long (panel) format with one row per
#'   subject-time combination.
#' @param id Character string giving the name of the subject identifier
#'   column. Default is \code{"id"}.
#' @param outcome Character string giving the name of the outcome column.
#'   Default is \code{"Y"}.
#' @param stat Character string specifying the summary statistic to apply
#'   to the outcome within each subject. One of \code{"min"}, \code{"max"},
#'   \code{"mean"}, or \code{"median"}. Default is \code{"mean"}.
#' @param formula A one-sided \code{\link[stats]{formula}} specifying the
#'   covariates for the linear model (e.g., \code{~ X + Z}). The outcome
#'   (as summarised) is used as the response. Default is \code{~ 1}
#'   (intercept only).
#' @param cov_stat Summary statistic to apply to time-varying covariates
#'   within each subject before fitting. Can be a single character string
#'   (applied to all covariates) or a named list mapping covariate names
#'   to statistic strings (e.g.,
#'   \code{list(X = "mean", W = "median")}). Covariates not named in
#'   the list default to \code{"mean"}. Default is \code{"mean"}.
#'
#' @return A numeric matrix of regression coefficients, identical to the
#'   \code{coefficients} component of
#'   \code{summary(\link[stats]{lm}(...))}, with columns \code{Estimate},
#'   \code{Std. Error}, \code{t value}, and \code{Pr(>|t|)}.
#'
#' @details
#' The function proceeds as follows:
#' \enumerate{
#'   \item The outcome is aggregated to one value per subject using
#'     \code{stat} (e.g., subject-level means).
#'   \item Each covariate in the formula is aggregated to one value per
#'     subject using \code{cov_stat}.
#'   \item A standard linear model is fitted on the subject-level
#'     summary data.
#' }
#'
#' This approach discards all within-subject temporal information and
#' does not account for informative visiting or observation processes.
#' It is included primarily as a naive benchmark for comparison with
#' methods such as \code{\link{GIVEHR}}, \code{\link{JMVL_IP}}, or
#' \code{\link{IIRR_weighting}}.
#'
#' @references
#' Fitzmaurice, G. M., Laird, N. M., and Ware, J. H. (2011).
#' \emph{Applied Longitudinal Analysis} (2nd ed.). Wiley. Chapter 4.
#'
#' @seealso
#' \code{\link{LMM}} for a mixed-model approach that uses all
#' observations;
#' \code{\link{GIVEHR}} for a joint model accounting for informative
#' visiting.
#'
#' @examples
#' \dontrun{
#' dat <- sim_data_gen(n = 500, seed = 42)
#'
#' # Mean outcome regressed on Z and X
#' SummaryStat(
#'   data    = dat$long_data,
#'   outcome = "Y",
#'   stat    = "mean",
#'   formula = ~ Z + X
#' )
#'
#' # Median outcome with per-covariate summarisation
#' SummaryStat(
#'   data     = dat$long_data,
#'   stat     = "median",
#'   formula  = ~ Z + X,
#'   cov_stat = list(X = "mean", Z = "mean")
#' )
#' }
#'
#' @importFrom stats aggregate as.formula lm
#' @export
SummaryStat <- function(data, id = "id", outcome = "Y", stat = "mean",
                        formula = ~ 1, cov_stat = "mean") {

  # Helper: get stat function
  get_stat_fun <- function(s) {
    switch(s,
           "min" = function(x) min(x, na.rm = TRUE),
           "max" = function(x) max(x, na.rm = TRUE),
           "mean" = function(x) mean(x, na.rm = TRUE),
           "median" = function(x) median(x, na.rm = TRUE),
           stop("stat must be 'min', 'max', 'mean', or 'median'"))
  }

  # Get covariates from formula
  covars <- all.vars(formula)

  # Summarize outcome by id
  outcome_agg <- aggregate(data[[outcome]], by = list(id = data[[id]]), FUN = get_stat_fun(stat))
  names(outcome_agg)[2] <- "outcome_summary"

  # Summarize covariates by id
  if (length(covars) > 0) {
    # Determine stat for each covariate
    cov_funs <- lapply(covars, function(v) {
      if (is.list(cov_stat)) {
        get_stat_fun(if (!is.null(cov_stat[[v]])) cov_stat[[v]] else "mean")
      } else {
        get_stat_fun(cov_stat)
      }
    })
    names(cov_funs) <- covars

    # Aggregate all covariates
    cov_agg_list <- lapply(covars, function(v) {
      aggregate(data[[v]], by = list(id = data[[id]]), FUN = cov_funs[[v]])$x
    })
    names(cov_agg_list) <- covars

    # Add to data frame
    outcome_agg <- cbind(outcome_agg, as.data.frame(cov_agg_list))
  }

  # Build formula and fit
  full_formula <- as.formula(paste("outcome_summary", paste(deparse(formula), collapse = "")))
  fit <- lm(full_formula, data = outcome_agg)

  return(summary(fit)$coefficients)
}

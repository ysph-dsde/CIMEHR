#' Summarize Longitudinal Data and Fit Linear Model
#'
#' @param data A data frame containing longitudinal data.
#' @param id Character string for subject ID column. Default is "id".
#' @param outcome Character string for outcome variable. Default is "Y".
#' @param stat Summary statistic for outcome: "min", "max", "mean", or "median". Default is "mean".
#' @param formula Formula for linear model (e.g., ~ X + Z).
#' @param cov_stat Summary statistic for time-varying covariates. Default is "mean".
#'        Can be a single value (applied to all) or a named list (e.g., list(X = "mean", W = "median")).
#'
#' @return Coefficient matrix from summary(lm(...))$coefficients.
#'
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

#dat = sim_ip_io_probit_re_lognormal(m = 1000, seed=1)
#dim(dat$long_data)

#fit = SummaryStat(data = data.frame(dat$long_data),
#                      id = "id",
#                      outcome = "Y",
#                      stat = "mean",
#                      formula = ~ Z + X)


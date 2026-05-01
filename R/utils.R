# ============================================================================
# Shared internal helpers
#
# This file consolidates utilities used by multiple package functions, to
# avoid duplicating identical definitions across files.
# ============================================================================


#' Null-coalescing operator
#'
#' Returns \code{y} if \code{x} is \code{NULL}, otherwise \code{x}.
#'
#' @param x,y Arbitrary objects.
#' @return \code{x} if non-\code{NULL}, else \code{y}.
#'
#' @keywords internal
#' @noRd
`%||%` <- function(x, y) if (is.null(x)) y else x


#' Validate Covariate Column Names Against a Data Frame
#'
#' Used by all CIMEHR fitting functions to give a clear error when a
#' user-supplied covariate name is not a column of \code{data}. The default
#' for every \code{covars_*} / \code{*_covs} argument across the package is
#' \code{NULL}, which is treated as "no covariates of this kind". This
#' helper standardises the validation and the error message so users get a
#' consistent hint across the package, including the list of available
#' columns.
#'
#' @param vars Character vector (or \code{NULL}) of column names supplied by
#'   the user.
#' @param data The \code{data.frame} the names should resolve in.
#' @param arg_name Name of the argument, used to make the error message
#'   informative.
#'
#' @return The input \code{vars} unchanged on success, with \code{NULL} and
#'   zero-length input both returned as \code{character(0)}.
#'
#' @keywords internal
#' @noRd
.check_covars <- function(vars, data, arg_name) {
  if (is.null(vars) || length(vars) == 0L) return(character(0))
  if (!is.character(vars)) {
    stop(sprintf("`%s` must be a character vector or NULL.", arg_name),
         call. = FALSE)
  }
  bad <- setdiff(vars, names(data))
  if (length(bad) > 0L) {
    stop(
      sprintf("`%s` references columns not found in `data`: %s.\n",
              arg_name, paste(sQuote(bad), collapse = ", ")),
      sprintf("Available columns: %s.\n",
              paste(sQuote(names(data)), collapse = ", ")),
      sprintf("Hint: pass the relevant columns explicitly, e.g. %s = c(\"Age\", \"Gender\").",
              arg_name),
      call. = FALSE
    )
  }
  vars
}


#' Apply an Informative Presence Adjustment Method
#'
#' Internal dispatcher that routes to the appropriate IP method function.
#' Used by both \code{Linear_increment_IP} and \code{Multiple_imputation_IP}.
#'
#' @param data A completed \code{data.frame} (after imputation).
#' @param ip_method Character string naming the method.
#' @param id,time,outcome,obs_indicator Column name strings.
#' @param visit_covs,outcome_covs,random_covs Character vectors of covariate
#'   names.
#' @param outcome_formula,random_formula Formula objects.
#' @param ip_args Named list of additional method-specific arguments.
#'
#' @return A list with elements \code{coef} and \code{se}.
#'
#' @keywords internal
#' @noRd
.apply_ip_method <- function(data, ip_method, id, time, outcome, obs_indicator,
                             visit_covs, outcome_covs, random_covs,
                             outcome_formula, random_formula,
                             ip_args = list()) {

  result <- switch(ip_method,
                   "Linear_mixed_model" = {
                     # Extract LMM-specific arguments with defaults
                     lmm_method <- ip_args$method %||% "REML"
                     lmm_optimizer <- ip_args$optimizer %||% "nlminb"
                     lmm_visit_adjust <- ip_args$visit_adjust  # NULL is valid default

                     coefs <- Linear_mixed_model(
                       data = data,
                       id = id,
                       fixed = outcome_formula,
                       random = random_formula,
                       obs_indicator = obs_indicator,
                       visit_adjust = lmm_visit_adjust,
                       method = lmm_method,
                       optimizer = lmm_optimizer
                     )
                     list(
                       coef = coefs$coefficients[, "Value"],
                       se = coefs$coefficients[, "Std.Error"]
                     )
                   },

                   "Joint_modeling_visiting_and_longitudinal_Liang" = {
                     fit <- Joint_modeling_visiting_and_longitudinal_Liang(
                       data = data,
                       id = id,
                       time = time,
                       outcome = outcome,
                       obs_indicator = obs_indicator,
                       visit_covs = visit_covs,
                       long_covs = outcome_covs,
                       random_covs = random_covs
                     )
                     list(
                       coef = fit$beta_hat,
                       se = NULL,  # Liang doesn't return SE directly
                       gamma_hat = fit$gamma_hat
                     )
                   },

                   "Inverse_intensity_rate_ratio_weighting" = {
                     # Extract IIRR-specific arguments with defaults
                     phi <- ip_args$phi %||% 0
                     balance_covs <- ip_args$balance_covs  # NULL uses visit_covs in function

                     fit <- Inverse_intensity_rate_ratio_balancing(
                       data = data,
                       id = id,
                       time = time,
                       outcome = outcome,
                       obs_indicator = obs_indicator,
                       visit_covs = visit_covs,
                       balance_covs = balance_covs,
                       formula = outcome_formula,
                       phi = phi
                     )
                     # Use MLE weights (standard IIRR)
                     list(
                       coef = fit$beta_mle,
                       se = NULL,
                       gamma_hat = fit$gamma_hat,
                       phi = phi
                     )
                   },

                   "Inverse_intensity_rate_ratio_balancing" = {
                     # Extract IIRR-specific arguments with defaults
                     phi <- ip_args$phi %||% 0
                     balance_covs <- ip_args$balance_covs  # NULL uses visit_covs in function

                     fit <- Inverse_intensity_rate_ratio_balancing(
                       data = data,
                       id = id,
                       time = time,
                       outcome = outcome,
                       obs_indicator = obs_indicator,
                       visit_covs = visit_covs,
                       balance_covs = balance_covs,
                       formula = outcome_formula,
                       phi = phi
                     )
                     # Use balancing weights
                     list(
                       coef = fit$beta_bal,
                       se = NULL,
                       gamma_hat = fit$gamma_hat,
                       phi = phi
                     )
                   },

                   "Pairwise_likelihood" = {
                     # Extract PairLik-specific arguments with defaults
                     pair_sample <- ip_args$pair_sample  # NULL uses all pairs
                     pair_family <- ip_args$family %||% "gaussian"

                     fit <- Pairwise_likelihood(
                       data = data[data[[obs_indicator]] == 1, ],
                       formula = outcome_formula,
                       id_col = id,
                       time_col = time,
                       y_col = outcome,
                       family = pair_family,
                       pair_sample = pair_sample
                     )
                     list(
                       coef = fit$coef,
                       se = if (!is.null(fit$vcov_naive)) sqrt(diag(fit$vcov_naive)) else NULL
                     )
                   },

                   stop("Unknown IP method: ", ip_method)
  )

  result
}# ============================================================================
# Shared internal helpers
#
# This file consolidates utilities used by multiple package functions, to
# avoid duplicating identical definitions across files.
# ============================================================================


#' Null-coalescing operator
#'
#' Returns \code{y} if \code{x} is \code{NULL}, otherwise \code{x}.
#'
#' @param x,y Arbitrary objects.
#' @return \code{x} if non-\code{NULL}, else \code{y}.
#'
#' @keywords internal
#' @noRd
`%||%` <- function(x, y) if (is.null(x)) y else x


#' Validate Covariate Column Names Against a Data Frame
#'
#' Used by all CIMEHR fitting functions to give a clear error when a
#' user-supplied covariate name is not a column of \code{data}. The default
#' for every \code{covars_*} / \code{*_covs} argument across the package is
#' \code{NULL}, which is treated as "no covariates of this kind". This
#' helper standardises the validation and the error message so users get a
#' consistent hint across the package, including the list of available
#' columns.
#'
#' @param vars Character vector (or \code{NULL}) of column names supplied by
#'   the user.
#' @param data The \code{data.frame} the names should resolve in.
#' @param arg_name Name of the argument, used to make the error message
#'   informative.
#'
#' @return The input \code{vars} unchanged on success, with \code{NULL} and
#'   zero-length input both returned as \code{character(0)}.
#'
#' @keywords internal
#' @noRd
.check_covars <- function(vars, data, arg_name) {
  if (is.null(vars) || length(vars) == 0L) return(character(0))
  if (!is.character(vars)) {
    stop(sprintf("`%s` must be a character vector or NULL.", arg_name),
         call. = FALSE)
  }
  bad <- setdiff(vars, names(data))
  if (length(bad) > 0L) {
    stop(
      sprintf("`%s` references columns not found in `data`: %s.\n",
              arg_name, paste(sQuote(bad), collapse = ", ")),
      sprintf("Available columns: %s.\n",
              paste(sQuote(names(data)), collapse = ", ")),
      sprintf("Hint: pass the relevant columns explicitly, e.g. %s = c(\"Age\", \"Gender\").",
              arg_name),
      call. = FALSE
    )
  }
  vars
}


#' Apply an Informative Presence Adjustment Method
#'
#' Internal dispatcher that routes to the appropriate IP method function.
#' Used by both \code{Linear_increment_IP} and \code{Multiple_imputation_IP}.
#'
#' @param data A completed \code{data.frame} (after imputation).
#' @param ip_method Character string naming the method.
#' @param id,time,outcome,obs_indicator Column name strings.
#' @param visit_covs,outcome_covs,random_covs Character vectors of covariate
#'   names.
#' @param outcome_formula,random_formula Formula objects.
#' @param ip_args Named list of additional method-specific arguments.
#'
#' @return A list with elements \code{coef} and \code{se}.
#'
#' @keywords internal
#' @noRd
.apply_ip_method <- function(data, ip_method, id, time, outcome, obs_indicator,
                             visit_covs, outcome_covs, random_covs,
                             outcome_formula, random_formula,
                             ip_args = list()) {

  result <- switch(ip_method,
                   "Linear_mixed_model" = {
                     # Extract LMM-specific arguments with defaults
                     lmm_method <- ip_args$method %||% "REML"
                     lmm_optimizer <- ip_args$optimizer %||% "nlminb"
                     lmm_visit_adjust <- ip_args$visit_adjust  # NULL is valid default

                     coefs <- Linear_mixed_model(
                       data = data,
                       id = id,
                       fixed = outcome_formula,
                       random = random_formula,
                       obs_indicator = obs_indicator,
                       visit_adjust = lmm_visit_adjust,
                       method = lmm_method,
                       optimizer = lmm_optimizer
                     )
                     list(
                       coef = coefs$coefficients[, "Value"],
                       se = coefs$coefficients[, "Std.Error"]
                     )
                   },

                   "Joint_modeling_visiting_and_longitudinal_Liang" = {
                     fit <- Joint_modeling_visiting_and_longitudinal_Liang(
                       data = data,
                       id = id,
                       time = time,
                       outcome = outcome,
                       obs_indicator = obs_indicator,
                       visit_covs = visit_covs,
                       long_covs = outcome_covs,
                       random_covs = random_covs
                     )
                     list(
                       coef = fit$beta_hat,
                       se = NULL,  # Liang doesn't return SE directly
                       gamma_hat = fit$gamma_hat
                     )
                   },

                   "Inverse_intensity_rate_ratio_weighting" = {
                     # Extract IIRR-specific arguments with defaults
                     phi <- ip_args$phi %||% 0
                     balance_covs <- ip_args$balance_covs  # NULL uses visit_covs in function

                     fit <- Inverse_intensity_rate_ratio_balancing(
                       data = data,
                       id = id,
                       time = time,
                       outcome = outcome,
                       obs_indicator = obs_indicator,
                       visit_covs = visit_covs,
                       balance_covs = balance_covs,
                       formula = outcome_formula,
                       phi = phi
                     )
                     # Use MLE weights (standard IIRR)
                     list(
                       coef = fit$beta_mle,
                       se = NULL,
                       gamma_hat = fit$gamma_hat,
                       phi = phi
                     )
                   },

                   "Inverse_intensity_rate_ratio_balancing" = {
                     # Extract IIRR-specific arguments with defaults
                     phi <- ip_args$phi %||% 0
                     balance_covs <- ip_args$balance_covs  # NULL uses visit_covs in function

                     fit <- Inverse_intensity_rate_ratio_balancing(
                       data = data,
                       id = id,
                       time = time,
                       outcome = outcome,
                       obs_indicator = obs_indicator,
                       visit_covs = visit_covs,
                       balance_covs = balance_covs,
                       formula = outcome_formula,
                       phi = phi
                     )
                     # Use balancing weights
                     list(
                       coef = fit$beta_bal,
                       se = NULL,
                       gamma_hat = fit$gamma_hat,
                       phi = phi
                     )
                   },

                   "Pairwise_likelihood" = {
                     # Extract PairLik-specific arguments with defaults
                     pair_sample <- ip_args$pair_sample  # NULL uses all pairs
                     pair_family <- ip_args$family %||% "gaussian"

                     fit <- Pairwise_likelihood(
                       data = data[data[[obs_indicator]] == 1, ],
                       formula = outcome_formula,
                       id_col = id,
                       time_col = time,
                       y_col = outcome,
                       family = pair_family,
                       pair_sample = pair_sample
                     )
                     list(
                       coef = fit$coef,
                       se = if (!is.null(fit$vcov_naive)) sqrt(diag(fit$vcov_naive)) else NULL
                     )
                   },

                   stop("Unknown IP method: ", ip_method)
  )

  result
}

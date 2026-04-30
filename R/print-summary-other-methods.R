# ============================================================================
# Print and summary methods for benchmark / comparison methods
#
# Methods covered:
#   - EHRJoint
#   - Joint_modeling_visiting_and_longitudinal_Liang
#   - Inverse_intensity_rate_ratio_weighting
#   - Inverse_intensity_rate_ratio_balancing
#   - Linear_mixed_model            (class "Linear_mixed_model_result")
#   - Summary_stat                  (class "Summary_stat_result")
#   - Pairwise_likelihood
#   - Linear_increment_IP           (class "Linear_increment_IP_result")
#   - Multiple_imputation_IP        (class "Multiple_imputation_IP_result")
#
# All methods share the helpers .CIMEHR_coef_table, .print_section, and
# .resolve_specify defined in print-summary-CIMEHR.R.
# ============================================================================


#' Build a Coefficient Table from an nlme Fixed-Effects Data Frame
#'
#' Appends Wald CI columns to the data frame returned by
#' \code{nlme::summary.lme()$tTable}.
#'
#' @param coef_df Data frame or numeric matrix with columns \code{Value}
#'   and \code{Std.Error}.
#' @param level Confidence level.
#'
#' @return The input data frame with additional columns \code{CI.low} and
#'   \code{CI.high}.
#'
#' @keywords internal
#' @noRd
.lmm_coef_table <- function(coef_df, level = 0.95) {
  coef_df <- as.data.frame(coef_df, check.names = FALSE)

  if (!all(c("Value", "Std.Error") %in% colnames(coef_df))) {
    stop("coef_df must contain columns 'Value' and 'Std.Error'.", call. = FALSE)
  }

  z <- stats::qnorm(1 - (1 - level) / 2)
  coef_df$CI.low  <- coef_df[["Value"]] - z * coef_df[["Std.Error"]]
  coef_df$CI.high <- coef_df[["Value"]] + z * coef_df[["Std.Error"]]
  coef_df
}


# ============================================================================
# EHRJoint
# ============================================================================

#' Print a Joint Estimation for Electronic Health Record (EHR) Longitudinal
#' Processes (EHRJoint) Fit
#'
#' Prints outcome fixed-effect estimates from the joint model.
#'
#' @param x An \code{EHRJoint} object returned by \code{\link{EHRJoint}}.
#' @param digits Number of significant digits. Default is 4.
#' @param level Confidence level for Wald intervals. Default is 0.95.
#' @param \dots Ignored.
#'
#' @return \code{x}, invisibly.
#'
#' @seealso \code{\link{summary.EHRJoint}}, \code{\link{EHRJoint}}
#'
#' @examples
#' dat <- sim_ehr_data[sim_ehr_data$id <= 80, ]
#' fit <- EHRJoint(dat)
#' print(fit)
#'
#' @export
print.EHRJoint <- function(x, digits = 4, level = 0.95, ...) {
  cat("\nEHRJoint -- Outcome Fixed Effect (beta)\n")
  cat(strrep("=", 50), "\n")

  tbl <- .CIMEHR_coef_table(x$beta_hat, se = NULL, level = level)
  print(round(tbl, digits))

  cat("\nNote: Standard errors are not available for EHRJoint.",
      "Consider bootstrap for inference.\n")
  cat("Use summary() for all process estimates, or\n")
  cat("summary(x, specify = \"visiting\" / \"observation\" / \"outcome\") for one layer.\n")
  invisible(x)
}

#' Summarise a Joint Estimation for Electronic Health Record (EHR)
#' Longitudinal Processes (EHRJoint) Fit
#'
#' Displays estimates for the requested layer(s): visiting, observation,
#' and/or outcome.
#'
#' @param object An \code{EHRJoint} object.
#' @param level Confidence level for Wald intervals. Default is 0.95.
#' @param digits Number of significant digits. Default is 4.
#' @param specify Character; which process layer to display. One of
#'   \code{"all"} (default), \code{"visiting"}, \code{"observation"}, or
#'   \code{"outcome"}.
#' @param \dots Ignored.
#'
#' @return \code{object}, invisibly.
#'
#' @seealso \code{\link{print.EHRJoint}}, \code{\link{EHRJoint}}
#'
#' @examples
#' dat <- sim_ehr_data[sim_ehr_data$id <= 80, ]
#' fit <- EHRJoint(dat)
#' summary(fit)
#' summary(fit, specify = "visiting")
#'
#' @export
summary.EHRJoint <- function(object, level = 0.95, digits = 4,
                             specify = c("all", "visiting", "observation", "outcome"),
                             ...) {
  specify <- match.arg(specify)
  layers  <- .resolve_specify(specify,
                              c("visiting", "observation", "outcome"))
  width   <- 60

  cat("\nEHRJoint Model -- Summary\n")
  cat(strrep("=", width), "\n")
  if (specify != "all")
    cat(sprintf("Layer requested: %s\n", specify))
  cat("\nNote: Standard errors are not available for this method.\n")
  cat("Consider bootstrap for inference.\n")

  if ("visiting" %in% layers) {
    .print_section("Visiting Process -- Fixed Effect (gamma)",
                   .CIMEHR_coef_table(object$gamma_hat, se = NULL, level),
                   digits, width)
    if (!is.null(object$sigma_sq_hat))
      cat(sprintf("\n  Frailty variance: %.4f\n", object$sigma_sq_hat))
  }

  if ("observation" %in% layers) {
    .print_section("Observation Process -- Fixed Effect (alpha)",
                   .CIMEHR_coef_table(object$alpha_hat, se = NULL, level),
                   digits, width)
  }

  if ("outcome" %in% layers) {
    .print_section("Longitudinal Outcome -- Fixed Effect (beta)",
                   .CIMEHR_coef_table(object$beta_hat, se = NULL, level),
                   digits, width)
    if (!is.null(object$theta_hat) && !all(is.na(object$theta_hat))) {
      .print_section("Longitudinal Outcome -- Random (Link) (theta)",
                     .CIMEHR_coef_table(object$theta_hat, se = NULL, level),
                     digits, width)
    }
  }

  cat("\n")
  invisible(object)
}


# ============================================================================
# Inverse_intensity_rate_ratio_weighting
# ============================================================================

#' Print an Inverse Intensity Rate Ratio (IIRR) Weighting Estimator Fit
#'
#' Prints outcome fixed-effect estimates from the inverse intensity rate
#' ratio weighting estimator.
#'
#' @param x An \code{Inverse_intensity_rate_ratio_weighting} object.
#' @param digits Number of significant digits. Default is 4.
#' @param level Confidence level for Wald intervals. Default is 0.95.
#' @param \dots Ignored.
#'
#' @return \code{x}, invisibly.
#'
#' @seealso \code{\link{summary.Inverse_intensity_rate_ratio_weighting}},
#'   \code{\link{Inverse_intensity_rate_ratio_weighting}}
#'
#' @examples
#' dat <- sim_ehr_data[sim_ehr_data$id <= 100, ]
#' fit <- Inverse_intensity_rate_ratio_weighting(dat, visit_covs = c("Z", "X"))
#' print(fit)
#'
#' @export
print.Inverse_intensity_rate_ratio_weighting <- function(x, digits = 4, level = 0.95, ...) {
  cat("\nIIRR Weighting -- Outcome Fixed Effect (beta)\n")
  cat(strrep("=", 50), "\n")

  tbl <- .CIMEHR_coef_table(x$beta_hat, se = NULL, level = level)
  print(round(tbl, digits))

  cat("\nNote: Standard errors are not available for",
      "Inverse_intensity_rate_ratio_weighting.\n")
  cat("Consider bootstrap for inference.\n")
  cat("Use summary() for all process estimates, or\n")
  cat("summary(x, specify = \"visiting\" / \"outcome\") for one layer.\n")
  invisible(x)
}

#' Summarise an Inverse Intensity Rate Ratio (IIRR) Weighting Estimator Fit
#'
#' Displays visiting-process and outcome estimates.
#'
#' @param object An \code{Inverse_intensity_rate_ratio_weighting} object.
#' @param level Confidence level for Wald intervals. Default is 0.95.
#' @param digits Number of significant digits. Default is 4.
#' @param specify Character; which process layer to display. One of
#'   \code{"all"} (default), \code{"visiting"}, or \code{"outcome"}.
#' @param \dots Ignored.
#'
#' @return \code{object}, invisibly.
#'
#' @seealso \code{\link{print.Inverse_intensity_rate_ratio_weighting}},
#'   \code{\link{Inverse_intensity_rate_ratio_weighting}}
#'
#' @examples
#' dat <- sim_ehr_data[sim_ehr_data$id <= 100, ]
#' fit <- Inverse_intensity_rate_ratio_weighting(dat, visit_covs = c("Z", "X"))
#' summary(fit)
#' summary(fit, specify = "visiting")
#'
#' @export
summary.Inverse_intensity_rate_ratio_weighting <- function(object, level = 0.95, digits = 4,
                                                           specify = c("all", "visiting", "outcome"),
                                                           ...) {
  specify <- match.arg(specify)
  layers  <- .resolve_specify(specify, c("visiting", "outcome"))
  width   <- 60

  cat("\nIIRR Weighting Estimator -- Summary\n")
  cat(strrep("=", width), "\n")
  if (specify != "all")
    cat(sprintf("Layer requested: %s\n", specify))
  cat("\nNote: Standard errors are not available for this method.\n")
  cat("Consider bootstrap for inference.\n")

  if ("visiting" %in% layers) {
    .print_section("Visiting Process -- Fixed Effect (gamma)",
                   .CIMEHR_coef_table(object$gamma_hat, se = NULL, level),
                   digits, width)
  }
  if ("outcome" %in% layers) {
    .print_section("Longitudinal Outcome -- Fixed Effect (beta)",
                   .CIMEHR_coef_table(object$beta_hat, se = NULL, level),
                   digits, width)
  }

  cat("\n")
  invisible(object)
}


# ============================================================================
# Inverse_intensity_rate_ratio_balancing
# ============================================================================

#' Print an Inverse Intensity Rate Ratio (IIRR) Estimator with Balancing
#' Weights Fit
#'
#' Prints outcome estimates from both balancing and maximum likelihood
#' estimation (MLE) weights.
#'
#' @param x An \code{Inverse_intensity_rate_ratio_balancing} object.
#' @param digits Number of significant digits. Default is 4.
#' @param level Confidence level for Wald intervals. Default is 0.95.
#' @param \dots Ignored.
#'
#' @return \code{x}, invisibly.
#'
#' @seealso \code{\link{summary.Inverse_intensity_rate_ratio_balancing}},
#'   \code{\link{Inverse_intensity_rate_ratio_balancing}}
#'
#' @examples
#' dat <- sim_ehr_data[sim_ehr_data$id <= 100, ]
#' fit <- Inverse_intensity_rate_ratio_balancing(dat, visit_covs = c("Z", "X"))
#' print(fit)
#'
#' @export
print.Inverse_intensity_rate_ratio_balancing <- function(x, digits = 4, level = 0.95, ...) {
  cat("\nIIRR Balancing -- Outcome Fixed Effect (beta)\n")
  cat(strrep("=", 50), "\n")

  cat("\nBalancing weights:\n")
  print(round(.CIMEHR_coef_table(x$beta_bal, se = NULL, level = level), digits))

  cat("\nMLE weights:\n")
  print(round(.CIMEHR_coef_table(x$beta_mle, se = NULL, level = level), digits))

  cat("\nNote: Standard errors are not available for",
      "Inverse_intensity_rate_ratio_balancing.\n")
  cat("Consider bootstrap for inference.\n")
  cat("Use summary() for all process estimates, or\n")
  cat("summary(x, specify = \"visiting\" / \"outcome\") for one layer.\n")
  invisible(x)
}

#' Summarise an Inverse Intensity Rate Ratio (IIRR) Estimator with
#' Balancing Weights Fit
#'
#' Displays visiting-process and outcome estimates under both balancing
#' and MLE weights, including the sensitivity parameter.
#'
#' @param object An \code{Inverse_intensity_rate_ratio_balancing} object.
#' @param level Confidence level for Wald intervals. Default is 0.95.
#' @param digits Number of significant digits. Default is 4.
#' @param specify Character; which process layer to display. One of
#'   \code{"all"} (default), \code{"visiting"}, or \code{"outcome"}.
#' @param \dots Ignored.
#'
#' @return \code{object}, invisibly.
#'
#' @seealso \code{\link{print.Inverse_intensity_rate_ratio_balancing}},
#'   \code{\link{Inverse_intensity_rate_ratio_balancing}}
#'
#' @examples
#' dat <- sim_ehr_data[sim_ehr_data$id <= 100, ]
#' fit <- Inverse_intensity_rate_ratio_balancing(dat, visit_covs = c("Z", "X"))
#' summary(fit)
#' summary(fit, specify = "outcome")
#'
#' @export
summary.Inverse_intensity_rate_ratio_balancing <- function(object, level = 0.95, digits = 4,
                                                           specify = c("all", "visiting", "outcome"),
                                                           ...) {
  specify <- match.arg(specify)
  layers  <- .resolve_specify(specify, c("visiting", "outcome"))
  width   <- 60

  cat("\nIIRR Balancing-Weights Estimator -- Summary\n")
  cat(strrep("=", width), "\n")
  if (specify != "all")
    cat(sprintf("Layer requested: %s\n", specify))
  cat("\nNote: Standard errors are not available for this method.\n")
  cat("Consider bootstrap for inference.\n")

  if (!is.null(object$phi))
    cat(sprintf("\n  Sensitivity parameter (phi): %.4f\n", object$phi))

  if ("visiting" %in% layers) {
    .print_section("Visiting Process -- Fixed Effect (gamma)",
                   .CIMEHR_coef_table(object$gamma_hat, se = NULL, level),
                   digits, width)
  }
  if ("outcome" %in% layers) {
    .print_section("Outcome -- Fixed Effect (beta) -- Balancing Weights",
                   .CIMEHR_coef_table(object$beta_bal, se = NULL, level),
                   digits, width)
    .print_section("Outcome -- Fixed Effect (beta) -- MLE Weights",
                   .CIMEHR_coef_table(object$beta_mle, se = NULL, level),
                   digits, width)
  }

  cat("\n")
  invisible(object)
}


# ============================================================================
# Joint_modeling_visiting_and_longitudinal_Liang
# ============================================================================

#' Print a Joint Model for Visiting and Longitudinal (JMVL) Processes
#' (Liang) Fit
#'
#' Prints outcome fixed-effect estimates from the Liang--Lu--Ying joint
#' estimating equations approach.
#'
#' @param x A \code{Joint_modeling_visiting_and_longitudinal_Liang} object.
#' @param digits Number of significant digits. Default is 4.
#' @param level Confidence level for Wald intervals. Default is 0.95.
#' @param \dots Ignored.
#'
#' @return \code{x}, invisibly.
#'
#' @seealso
#'   \code{\link{summary.Joint_modeling_visiting_and_longitudinal_Liang}},
#'   \code{\link{Joint_modeling_visiting_and_longitudinal_Liang}}
#'
#' @examples
#' dat <- sim_ehr_data[sim_ehr_data$id <= 80, ]
#' fit <- Joint_modeling_visiting_and_longitudinal_Liang(dat)
#' print(fit)
#'
#' @export
print.Joint_modeling_visiting_and_longitudinal_Liang <- function(x, digits = 4, level = 0.95, ...) {
  cat("\nJMVL (Liang) -- Outcome Fixed Effect (beta)\n")
  cat(strrep("=", 50), "\n")

  tbl <- .CIMEHR_coef_table(x$beta_hat, se = NULL, level = level)
  print(round(tbl, digits))

  cat("\nNote: Standard errors are not available for",
      "Joint_modeling_visiting_and_longitudinal_Liang.\n")
  cat("Consider bootstrap for inference.\n")
  cat("Use summary() for all process estimates, or\n")
  cat("summary(x, specify = \"visiting\" / \"outcome\") for one layer.\n")
  invisible(x)
}

#' Summarise a Joint Model for Visiting and Longitudinal (JMVL) Processes
#' (Liang) Fit
#'
#' Displays visiting-process and outcome estimates.
#'
#' @param object A \code{Joint_modeling_visiting_and_longitudinal_Liang}
#'   object.
#' @param level Confidence level for Wald intervals. Default is 0.95.
#' @param digits Number of significant digits. Default is 4.
#' @param specify Character; which process layer to display. One of
#'   \code{"all"} (default), \code{"visiting"}, or \code{"outcome"}.
#' @param \dots Ignored.
#'
#' @return \code{object}, invisibly.
#'
#' @seealso \code{\link{print.Joint_modeling_visiting_and_longitudinal_Liang}},
#'   \code{\link{Joint_modeling_visiting_and_longitudinal_Liang}}
#'
#' @examples
#' dat <- sim_ehr_data[sim_ehr_data$id <= 80, ]
#' fit <- Joint_modeling_visiting_and_longitudinal_Liang(dat)
#' summary(fit)
#' summary(fit, specify = "visiting")
#'
#' @export
summary.Joint_modeling_visiting_and_longitudinal_Liang <- function(object, level = 0.95, digits = 4,
                                                                   specify = c("all", "visiting", "outcome"),
                                                                   ...) {
  specify <- match.arg(specify)
  layers  <- .resolve_specify(specify, c("visiting", "outcome"))
  width   <- 60

  cat("\nJMVL (Liang) Joint Model -- Summary\n")
  cat(strrep("=", width), "\n")
  if (specify != "all")
    cat(sprintf("Layer requested: %s\n", specify))
  cat("\nNote: Standard errors are not available for this method.\n")
  cat("Consider bootstrap for inference.\n")

  if ("visiting" %in% layers) {
    .print_section("Visiting Process -- Fixed Effect (gamma)",
                   .CIMEHR_coef_table(object$gamma_hat, se = NULL, level),
                   digits, width)
    if (!is.null(object$sigma_sq_hat))
      cat(sprintf("\n  Frailty variance: %.4f\n", object$sigma_sq_hat))
  }
  if ("outcome" %in% layers) {
    .print_section("Longitudinal Outcome -- Fixed Effect (beta)",
                   .CIMEHR_coef_table(object$beta_hat, se = NULL, level),
                   digits, width)
  }

  cat("\n")
  invisible(object)
}


# ============================================================================
# Linear_mixed_model_result
# ============================================================================

#' Print a Linear Mixed Model for Longitudinal Data Result
#'
#' Prints the fixed-effects table from a linear mixed model fit.
#'
#' @param x A \code{Linear_mixed_model_result} object returned by
#'   \code{\link{Linear_mixed_model}}.
#' @param digits Number of significant digits. Default is 4.
#' @param level Confidence level for Wald intervals. Default is 0.95.
#' @param \dots Ignored.
#'
#' @return \code{x}, invisibly.
#'
#' @seealso \code{\link{summary.Linear_mixed_model_result}},
#'   \code{\link{Linear_mixed_model}}
#'
#' @examples
#' dat <- sim_ehr_data[sim_ehr_data$id <= 100, ]
#' fit <- Linear_mixed_model(dat, fixed = Y ~ Z + X + time,
#'            random = ~ 1 | id, obs_indicator = "R")
#' print(fit)
#'
#' @export
print.Linear_mixed_model_result <- function(x, digits = 4, level = 0.95, ...) {
  cat("\nLinear Mixed Model -- Fixed Effect\n")
  cat(strrep("=", 50), "\n")

  tbl <- .lmm_coef_table(x$coefficients, level)
  print(round(tbl, digits))
  invisible(x)
}

#' Summarise a Linear Mixed Model for Longitudinal Data Result
#'
#' Displays the fixed-effects coefficient table.
#'
#' @param object A \code{Linear_mixed_model_result} object.
#' @param level Confidence level for Wald intervals. Default is 0.95.
#' @param digits Number of significant digits. Default is 4.
#' @param specify Character; only \code{"all"} or \code{"outcome"} are
#'   meaningful for this single-layer fit. Default is \code{"all"}.
#' @param \dots Ignored.
#'
#' @return \code{object}, invisibly.
#'
#' @seealso \code{\link{print.Linear_mixed_model_result}},
#'   \code{\link{Linear_mixed_model}}
#'
#' @examples
#' dat <- sim_ehr_data[sim_ehr_data$id <= 100, ]
#' fit <- Linear_mixed_model(dat, fixed = Y ~ Z + X + time,
#'            random = ~ 1 | id, obs_indicator = "R")
#' summary(fit)
#'
#' @export
summary.Linear_mixed_model_result <- function(object, level = 0.95, digits = 4,
                                              specify = c("all", "outcome"),
                                              ...) {
  specify <- match.arg(specify)
  width   <- 60

  cat("\nLinear Mixed Model -- Summary\n")
  cat(strrep("=", width), "\n")
  if (specify != "all")
    cat(sprintf("Layer requested: %s\n", specify))

  .print_section("Fixed Effect",
                 .lmm_coef_table(object$coefficients, level), digits, width)

  cat("\n")
  invisible(object)
}


# ============================================================================
# Summary_stat_result
# ============================================================================

#' Print a Summary-Statistic Regression for Longitudinal Data Result
#'
#' Prints coefficients from the summary-statistic regression.
#'
#' @param x A \code{Summary_stat_result} object returned by
#'   \code{\link{Summary_stat}}.
#' @param digits Number of significant digits. Default is 4.
#' @param \dots Ignored.
#'
#' @return \code{x}, invisibly.
#'
#' @seealso \code{\link{summary.Summary_stat_result}},
#'   \code{\link{Summary_stat}}
#'
#' @examples
#' dat <- sim_ehr_data[sim_ehr_data$id <= 100, ]
#' fit <- Summary_stat(dat, formula = ~ Z + X)
#' print(fit)
#'
#' @export
print.Summary_stat_result <- function(x, digits = 4, ...) {
  cat("\nSummary-Statistic Regression -- Coefficients (Fixed Effect)\n")
  cat(strrep("=", 50), "\n")

  print(round(x$coefficients, digits))
  invisible(x)
}

#' Summarise a Summary-Statistic Regression for Longitudinal Data Result
#'
#' Displays the regression coefficient table.
#'
#' @param object A \code{Summary_stat_result} object.
#' @param digits Number of significant digits. Default is 4.
#' @param specify Character; only \code{"all"} or \code{"outcome"} are
#'   meaningful for this single-layer fit. Default is \code{"all"}.
#' @param \dots Ignored.
#'
#' @return \code{object}, invisibly.
#'
#' @seealso \code{\link{print.Summary_stat_result}},
#'   \code{\link{Summary_stat}}
#'
#' @examples
#' dat <- sim_ehr_data[sim_ehr_data$id <= 100, ]
#' fit <- Summary_stat(dat, formula = ~ Z + X)
#' summary(fit)
#'
#' @export
summary.Summary_stat_result <- function(object, digits = 4,
                                        specify = c("all", "outcome"),
                                        ...) {
  specify <- match.arg(specify)
  width   <- 60

  cat("\nSummary-Statistic Regression -- Summary\n")
  cat(strrep("=", width), "\n")
  if (specify != "all")
    cat(sprintf("Layer requested: %s\n", specify))

  cat("\nFixed Effect (Coefficients):\n")
  cat(strrep("-", width), "\n")
  print(round(object$coefficients, digits))

  cat("\n")
  invisible(object)
}


# ============================================================================
# Pairwise_likelihood
# ============================================================================

#' Print a Pairwise Composite Likelihood (PCL) Generalized Linear Model
#' (Pairwise_likelihood) Fit
#'
#' Prints coefficients with naive standard errors (from the inverse
#' Hessian) when available.
#'
#' @param x A \code{Pairwise_likelihood} object returned by
#'   \code{\link{Pairwise_likelihood}}.
#' @param digits Number of significant digits. Default is 4.
#' @param level Confidence level for Wald intervals. Default is 0.95.
#' @param \dots Ignored.
#'
#' @return \code{x}, invisibly.
#'
#' @seealso \code{\link{summary.Pairwise_likelihood}},
#'   \code{\link{Pairwise_likelihood}}
#'
#' @examples
#' dat <- sim_ehr_data[sim_ehr_data$id <= 100 &
#'                     sim_ehr_data$R == 1 & !is.na(sim_ehr_data$R), ]
#' fit <- Pairwise_likelihood(dat, formula = Y ~ Z + X + time,
#'                            pair_sample = 200)
#' print(fit)
#'
#' @export
print.Pairwise_likelihood <- function(x, digits = 4, level = 0.95, ...) {
  cat("\nPairwise Likelihood -- Fixed Effect (Coefficients)\n")
  cat(strrep("=", 50), "\n")

  se <- NULL
  if (!is.null(x$vcov_naive))
    se <- sqrt(pmax(diag(x$vcov_naive), 0))

  tbl <- .CIMEHR_coef_table(x$coef, se, level = level)
  print(round(tbl, digits))

  cat(sprintf("\nConverged: %s | Log-likelihood: %.2f\n",
              x$converged, x$loglik))
  invisible(x)
}

#' Summarise a Pairwise Composite Likelihood (PCL) Generalized Linear
#' Model (Pairwise_likelihood) Fit
#'
#' Displays model dimensions, convergence status, and the coefficient
#' table with standard errors when available.
#'
#' @param object A \code{Pairwise_likelihood} object.
#' @param level Confidence level for Wald intervals. Default is 0.95.
#' @param digits Number of significant digits. Default is 4.
#' @param specify Character; only \code{"all"} or \code{"outcome"} are
#'   meaningful for this single-layer fit. Default is \code{"all"}.
#' @param \dots Ignored.
#'
#' @return \code{object}, invisibly.
#'
#' @seealso \code{\link{print.Pairwise_likelihood}},
#'   \code{\link{Pairwise_likelihood}}
#'
#' @examples
#' dat <- sim_ehr_data[sim_ehr_data$id <= 100 &
#'                     sim_ehr_data$R == 1 & !is.na(sim_ehr_data$R), ]
#' fit <- Pairwise_likelihood(dat, formula = Y ~ Z + X + time,
#'                            pair_sample = 200)
#' summary(fit)
#'
#' @export
summary.Pairwise_likelihood <- function(object, level = 0.95, digits = 4,
                                        specify = c("all", "outcome"),
                                        ...) {
  specify <- match.arg(specify)
  width   <- 60
  cat("\nPairwise Likelihood -- Summary\n")
  cat(strrep("=", width), "\n")
  if (specify != "all")
    cat(sprintf("Layer requested: %s\n", specify))

  cat(sprintf("\n  Observations: %d | Clusters: %d | Pairs: %d\n",
              object$n, object$G, object$k_pairs))
  cat(sprintf("  Converged: %s | Log-likelihood: %.2f\n",
              object$converged, object$loglik))

  se <- NULL
  if (!is.null(object$vcov_naive))
    se <- sqrt(pmax(diag(object$vcov_naive), 0))

  .print_section("Fixed Effect (Coefficients)",
                 .CIMEHR_coef_table(object$coef, se, level), digits, width)

  cat("\n")
  invisible(object)
}


# ============================================================================
# Linear_increment_IP_result
# ============================================================================

#' Print a Linear Increment (LI) Imputation with Informative Presence (IP)
#' Adjustment Fit
#'
#' Prints outcome coefficients from the linear increment imputation with
#' informative presence adjustment. Standard errors and confidence
#' intervals are not reported; see the bootstrap note in the output.
#'
#' @param x A \code{Linear_increment_IP_result} object returned by
#'   \code{\link{Linear_increment_IP}}.
#' @param digits Number of significant digits. Default is 4.
#' @param \dots Ignored.
#'
#' @return \code{x}, invisibly.
#'
#' @seealso \code{\link{summary.Linear_increment_IP_result}},
#'   \code{\link{Linear_increment_IP}}
#'
#' @examples
#' dat <- sim_ehr_data[sim_ehr_data$id <= 80,
#'                     c("id", "Z", "X", "time", "R", "Y", "C")]
#' fit <- Linear_increment_IP(data = dat, ip_method = "Linear_mixed_model",
#'              visit_covs = c("Z", "X"))
#' print(fit)
#'
#' @export
print.Linear_increment_IP_result <- function(x, digits = 4, ...) {
  cat("\nLinear Increment + IP -- Outcome Fixed Effect (Coefficients)\n")
  cat(strrep("=", 50), "\n")

  tbl <- .CIMEHR_coef_table(x$coefficients, se = NULL)
  print(round(tbl, digits))

  cat("\nNote: Standard errors and confidence intervals are not reported.\n")
  cat("Valid inference requires bootstrapping the full two-stage procedure\n")
  cat("(LI imputation + IP analysis) over resampled subjects, as SE from\n")
  cat("the IP model alone does not account for imputation uncertainty.\n")
  cat("We suggest users write their own bootstrap method.\n")

  cat(sprintf("\nIP method: %s | Values imputed: %d\n",
              x$ip_method, x$n_imputed))
  invisible(x)
}

#' Summarise a Linear Increment (LI) Imputation with Informative Presence
#' (IP) Adjustment Fit
#'
#' Displays imputation model details and outcome coefficients. Standard
#' errors and confidence intervals are not reported; see the bootstrap
#' note in the output.
#'
#' @param object A \code{Linear_increment_IP_result} object.
#' @param digits Number of significant digits. Default is 4.
#' @param specify Character; only \code{"all"} or \code{"outcome"} are
#'   meaningful for this two-stage fit (only the outcome stage is
#'   reported). Default is \code{"all"}.
#' @param \dots Ignored.
#'
#' @return \code{object}, invisibly.
#'
#' @seealso \code{\link{print.Linear_increment_IP_result}},
#'   \code{\link{Linear_increment_IP}}
#'
#' @examples
#' dat <- sim_ehr_data[sim_ehr_data$id <= 80,
#'                     c("id", "Z", "X", "time", "R", "Y", "C")]
#' fit <- Linear_increment_IP(data = dat, ip_method = "Linear_mixed_model",
#'              visit_covs = c("Z", "X"))
#' summary(fit)
#'
#' @export
summary.Linear_increment_IP_result <- function(object, digits = 4,
                                               specify = c("all", "outcome"),
                                               ...) {
  specify <- match.arg(specify)
  width   <- 60
  cat("\nLinear Increment Imputation + IP Adjustment -- Summary\n")
  cat(strrep("=", width), "\n")
  if (specify != "all")
    cat(sprintf("Layer requested: %s\n", specify))

  cat(sprintf("\n  LI formula: %s\n", deparse(object$li_model$formula)))
  cat(sprintf("  Covariance: %s | Values imputed: %d\n",
              object$li_model$covariance, object$n_imputed))
  cat(sprintf("  IP method: %s\n", object$ip_method))

  cat("\n  Note: Standard errors and confidence intervals are not reported.\n")
  cat("  Valid inference requires bootstrapping the full two-stage procedure\n")
  cat("  (LI imputation + IP analysis) over resampled subjects, as SE from\n")
  cat("  the IP model alone does not account for imputation uncertainty.\n")
  cat("  We suggest users write their own bootstrap method.\n")

  .print_section("Outcome -- Fixed Effect (Coefficients)",
                 .CIMEHR_coef_table(object$coefficients, se = NULL),
                 digits, width)

  cat("\n")
  invisible(object)
}


# ============================================================================
# Multiple_imputation_IP_result
# ============================================================================

#' Print a Multiple Imputation (MI) with Informative Presence (IP)
#' Adjustment Fit
#'
#' Prints pooled outcome coefficients from the multiple imputation with
#' informative presence adjustment. Standard errors and confidence
#' intervals are not reported; see the bootstrap note in the output.
#'
#' @param x A \code{Multiple_imputation_IP_result} object returned by
#'   \code{\link{Multiple_imputation_IP}}.
#' @param digits Number of significant digits. Default is 4.
#' @param \dots Ignored.
#'
#' @return \code{x}, invisibly.
#'
#' @seealso \code{\link{summary.Multiple_imputation_IP_result}},
#'   \code{\link{Multiple_imputation_IP}}
#'
#' @examples
#' dat <- sim_ehr_data[sim_ehr_data$id <= 80,
#'                     c("id", "Z", "X", "time", "R", "Y", "C")]
#' fit <- Multiple_imputation_IP(data = dat, ip_method = "Linear_mixed_model",
#'              visit_covs = c("Z", "X"), m = 2, maxit = 5)
#' print(fit)
#'
#' @export
print.Multiple_imputation_IP_result <- function(x, digits = 4, ...) {
  cat("\nMultiple Imputation + IP -- Outcome Fixed Effect (Coefficients)\n")
  cat(strrep("=", 50), "\n")

  tbl <- .CIMEHR_coef_table(x$coefficients, se = NULL)
  print(round(tbl, digits))

  cat("\nNote: Standard errors and confidence intervals are not reported.\n")
  cat("Rubin's rules require per-imputation SE from the IP method, which\n")
  cat("is not available for all IP methods; and even where available, the\n")
  cat("combined SE does not account for uncertainty in the imputation model\n")
  cat("itself. For valid inference, bootstrap the full procedure\n")
  cat("(imputation + IP analysis) over resampled subjects.\n")
  cat("We suggest users write their own bootstrap method.\n")

  cat(sprintf("\nIP method: %s | Imputations: %d/%d\n",
              x$ip_method, x$m, x$m_requested))
  invisible(x)
}

#' Summarise a Multiple Imputation (MI) with Informative Presence (IP)
#' Adjustment Fit
#'
#' Displays imputation settings and pooled outcome coefficients. Standard
#' errors and confidence intervals are not reported; see the bootstrap
#' note in the output.
#'
#' @param object A \code{Multiple_imputation_IP_result} object.
#' @param digits Number of significant digits. Default is 4.
#' @param specify Character; only \code{"all"} or \code{"outcome"} are
#'   meaningful for this two-stage fit (only the pooled outcome stage is
#'   reported). Default is \code{"all"}.
#' @param \dots Ignored.
#'
#' @return \code{object}, invisibly.
#'
#' @seealso \code{\link{print.Multiple_imputation_IP_result}},
#'   \code{\link{Multiple_imputation_IP}}
#'
#' @examples
#' dat <- sim_ehr_data[sim_ehr_data$id <= 80,
#'                     c("id", "Z", "X", "time", "R", "Y", "C")]
#' fit <- Multiple_imputation_IP(data = dat, ip_method = "Linear_mixed_model",
#'              visit_covs = c("Z", "X"), m = 2, maxit = 5)
#' summary(fit)
#'
#' @export
summary.Multiple_imputation_IP_result <- function(object, digits = 4,
                                                  specify = c("all", "outcome"),
                                                  ...) {
  specify <- match.arg(specify)
  width   <- 60
  cat("\nMultiple Imputation + IP Adjustment -- Summary\n")
  cat(strrep("=", width), "\n")
  if (specify != "all")
    cat(sprintf("Layer requested: %s\n", specify))

  cat(sprintf("\n  Imputation method: %s\n", object$impute_method))
  cat(sprintf("  Imputations used: %d / %d requested\n",
              object$m, object$m_requested))
  cat(sprintf("  IP method: %s\n", object$ip_method))

  cat("\n  Note: Standard errors and confidence intervals are not reported.\n")
  cat("  Rubin's rules require per-imputation SE from the IP method, which\n")
  cat("  is not available for all IP methods; and even where available, the\n")
  cat("  combined SE does not account for uncertainty in the imputation model\n")
  cat("  itself. For valid inference, bootstrap the full procedure\n")
  cat("  (imputation + IP analysis) over resampled subjects.\n")
  cat("  We suggest users write their own bootstrap method.\n")

  .print_section("Outcome -- Pooled Fixed Effect (Coefficients)",
                 .CIMEHR_coef_table(object$coefficients, se = NULL),
                 digits, width)

  cat("\n")
  invisible(object)
}

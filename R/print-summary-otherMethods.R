# ============================================================================
# Print and summary methods for benchmark / comparison methods
# ============================================================================


# ============================================================================
# Internal helpers
# ============================================================================

#' Build a Coefficient Table from an nlme Fixed-Effects Data Frame
#'
#' Appends Wald CI columns to the data frame returned by
#' \code{nlme::intervals()} or \code{nlme::summary.lme()$tTable}.
#'
#' @param coef_df Data frame or numeric matrix with columns \code{Value} and \code{Std.Error}.
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

#' Print an Joint Estimation for Electronic Health Record (EHR) Longitudinal
#' Processes (EHRJoint) Fit
#'
#' Prints outcome fixed-effect estimates from the joint model without a
#' separate observation process.
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
  cat("\nEHRJoint -- Outcome Fixed Effects (beta)\n")
  cat(strrep("=", 50), "\n")

  tbl <- .givehr_coef_table(x$beta_hat, se = NULL, level = level)
  print(round(tbl, digits))

  cat("\nNote: Standard errors are not available for EHRJoint.",
      "Consider bootstrap for inference.\n")
  cat("Use summary() for all process estimates.\n")
  invisible(x)
}

#' Summarise an Joint Estimation for Electronic Health Record (EHR) Longitudinal
#' Processes (EHRJoint) Fit
#'
#' Displays estimates for the visiting, observation, and outcome processes.
#'
#' @param object An \code{EHRJoint} object.
#' @param level Confidence level for Wald intervals. Default is 0.95.
#' @param digits Number of significant digits. Default is 4.
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
#'
#' @export
summary.EHRJoint <- function(object, level = 0.95, digits = 4, ...) {
  width <- 60
  cat("\nEHRJoint Model -- Summary\n")
  cat(strrep("=", width), "\n")
  cat("\nNote: Standard errors are not available for this method.\n")
  cat("Consider bootstrap for inference.\n")

  .print_section("Visiting Process (gamma)",
                 .givehr_coef_table(object$gamma_hat, se = NULL, level), digits, width)

  if (!is.null(object$sigma_sq_hat))
    cat(sprintf("\n  Frailty variance: %.4f\n", object$sigma_sq_hat))

  .print_section("Observation Process (alpha)",
                 .givehr_coef_table(object$alpha_hat, se = NULL, level), digits, width)

  .print_section("Longitudinal Outcome -- Fixed Effects (beta)",
                 .givehr_coef_table(object$beta_hat, se = NULL, level), digits, width)

  if (!is.null(object$theta_hat) && !all(is.na(object$theta_hat))) {
    .print_section("Longitudinal Outcome -- Latent Link (theta)",
                   .givehr_coef_table(object$theta_hat, se = NULL, level), digits, width)
  }

  cat("\n")
  invisible(object)
}


# ============================================================================
# IIRR_weighting
# ============================================================================

#' Print an Inverse Intensity Rate Ratio (IIRR) Weighting Estimator (IIRR_weighting) Fit
#'
#' Prints outcome fixed-effect estimates from the inverse intensity rate
#' ratio weighting estimator.
#'
#' @param x An \code{IIRR_weighting} object.
#' @param digits Number of significant digits. Default is 4.
#' @param level Confidence level for Wald intervals. Default is 0.95.
#' @param \dots Ignored.
#'
#' @return \code{x}, invisibly.
#'
#' @seealso \code{\link{summary.IIRR_weighting}},
#'   \code{\link{IIRR_weighting}}
#'
#' @examples
#' dat <- sim_ehr_data[sim_ehr_data$id <= 100, ]
#' fit <- IIRR_weighting(dat, visit_covs = c("Z", "X"))
#' print(fit)
#'
#' @export
print.IIRR_weighting <- function(x, digits = 4, level = 0.95, ...) {
  cat("\nIIRR Weighting -- Outcome Fixed Effects (beta)\n")
  cat(strrep("=", 50), "\n")

  tbl <- .givehr_coef_table(x$beta_hat, se = NULL, level = level)
  print(round(tbl, digits))

  cat("\nNote: Standard errors are not available for IIRR_weighting.",
      "Consider bootstrap for inference.\n")
  cat("Use summary() for all process estimates.\n")
  invisible(x)
}

#' Summarise an Inverse Intensity Rate Ratio (IIRR) Weighting Estimator (IIRR_weighting) Fit
#'
#' Displays visiting-process and outcome estimates.
#'
#' @param object An \code{IIRR_weighting} object.
#' @param level Confidence level for Wald intervals. Default is 0.95.
#' @param digits Number of significant digits. Default is 4.
#' @param \dots Ignored.
#'
#' @return \code{object}, invisibly.
#'
#' @seealso \code{\link{print.IIRR_weighting}},
#'   \code{\link{IIRR_weighting}}
#'
#' @examples
#' dat <- sim_ehr_data[sim_ehr_data$id <= 100, ]
#' fit <- IIRR_weighting(dat, visit_covs = c("Z", "X"))
#' summary(fit)
#'
#' @export
summary.IIRR_weighting <- function(object, level = 0.95, digits = 4, ...) {
  width <- 60
  cat("\nIIRR Weighting Estimator -- Summary\n")
  cat(strrep("=", width), "\n")
  cat("\nNote: Standard errors are not available for this method.\n")
  cat("Consider bootstrap for inference.\n")

  .print_section("Visiting Process (gamma)",
                 .givehr_coef_table(object$gamma_hat, se = NULL, level), digits, width)

  .print_section("Longitudinal Outcome -- Fixed Effects (beta)",
                 .givehr_coef_table(object$beta_hat, se = NULL, level), digits, width)

  cat("\n")
  invisible(object)
}


# ============================================================================
# IIRR_balancing
# ============================================================================

#' Print an Inverse Intensity Rate Ratio (IIRR) Estimator with Balancing Weights
#' (IIRR_balancing) Fit
#'
#' Prints outcome estimates from both balancing and maximum likelihood estimation (MLE) weights.
#'
#' @param x An \code{IIRR_balancing} object.
#' @param digits Number of significant digits. Default is 4.
#' @param level Confidence level for Wald intervals. Default is 0.95.
#' @param \dots Ignored.
#'
#' @return \code{x}, invisibly.
#'
#' @seealso \code{\link{summary.IIRR_balancing}},
#'   \code{\link{IIRR_balancing}}
#'
#' @examples
#' dat <- sim_ehr_data[sim_ehr_data$id <= 100, ]
#' fit <- IIRR_balancing(dat, visit_covs = c("Z", "X"))
#' print(fit)
#'
#' @export
print.IIRR_balancing <- function(x, digits = 4, level = 0.95, ...) {
  cat("\nIIRR Balancing -- Outcome Fixed Effects (beta)\n")
  cat(strrep("=", 50), "\n")

  cat("\nBalancing weights:\n")
  print(round(.givehr_coef_table(x$beta_bal, se = NULL, level = level), digits))

  cat("\nMLE weights:\n")
  print(round(.givehr_coef_table(x$beta_mle, se = NULL, level = level), digits))

  cat("\nNote: Standard errors are not available for IIRR_balancing.",
      "Consider bootstrap for inference.\n")
  cat("Use summary() for all process estimates.\n")
  invisible(x)
}

#' Summarise an Inverse Intensity Rate Ratio (IIRR) Estimator with Balancing Weights
#' (IIRR_balancing) Fit
#'
#' Displays visiting-process and outcome estimates under both balancing
#' and MLE weights, including the sensitivity parameter.
#'
#' @param object An \code{IIRR_balancing} object.
#' @param level Confidence level for Wald intervals. Default is 0.95.
#' @param digits Number of significant digits. Default is 4.
#' @param \dots Ignored.
#'
#' @return \code{object}, invisibly.
#'
#' @seealso \code{\link{print.IIRR_balancing}},
#'   \code{\link{IIRR_balancing}}
#'
#' @examples
#' dat <- sim_ehr_data[sim_ehr_data$id <= 100, ]
#' fit <- IIRR_balancing(dat, visit_covs = c("Z", "X"))
#' summary(fit)
#'
#' @export
summary.IIRR_balancing <- function(object, level = 0.95, digits = 4, ...) {
  width <- 60
  cat("\nIIRR Balancing-Weights Estimator -- Summary\n")
  cat(strrep("=", width), "\n")
  cat("\nNote: Standard errors are not available for this method.\n")
  cat("Consider bootstrap for inference.\n")

  if (!is.null(object$phi))
    cat(sprintf("\n  Sensitivity parameter (phi): %.4f\n", object$phi))

  .print_section("Visiting Process (gamma)",
                 .givehr_coef_table(object$gamma_hat, se = NULL, level), digits, width)

  .print_section("Outcome (beta) -- Balancing Weights",
                 .givehr_coef_table(object$beta_bal, se = NULL, level), digits, width)

  .print_section("Outcome (beta) -- MLE Weights",
                 .givehr_coef_table(object$beta_mle, se = NULL, level), digits, width)

  cat("\n")
  invisible(object)
}


# ============================================================================
# JMVL_Liang
# ============================================================================

#' Print a Joint Model for Visiting and Longitudinal (JMVL) Processes with
#' Informative Presence (IP) Adjustment (JMVL_Liang) Fit
#'
#' Prints outcome fixed-effect estimates from the Liang--Lu--Ying joint
#' estimating equations approach.
#'
#' @param x A \code{JMVL_Liang} object.
#' @param digits Number of significant digits. Default is 4.
#' @param level Confidence level for Wald intervals. Default is 0.95.
#' @param \dots Ignored.
#'
#' @return \code{x}, invisibly.
#'
#' @seealso \code{\link{summary.JMVL_Liang}}, \code{\link{JMVL_Liang}}
#'
#' @examples
#' dat <- sim_ehr_data[sim_ehr_data$id <= 80, ]
#' fit <- JMVL_Liang(dat)
#' print(fit)
#'
#' @export
print.JMVL_Liang <- function(x, digits = 4, level = 0.95, ...) {
  cat("\nJMVL_Liang -- Outcome Fixed Effects (beta)\n")
  cat(strrep("=", 50), "\n")

  tbl <- .givehr_coef_table(x$beta_hat, se = NULL, level = level)
  print(round(tbl, digits))

  cat("\nNote: Standard errors are not available for JMVL_Liang.",
      "Consider bootstrap for inference.\n")
  cat("Use summary() for all process estimates.\n")
  invisible(x)
}

#' Summarise a Joint Model for Visiting and Longitudinal (JMVL) Processes with
#' Informative Presence (IP) Adjustment (JMVL_Liang) Fit
#'
#' Displays visiting-process and outcome estimates.
#'
#' @param object A \code{JMVL_Liang} object.
#' @param level Confidence level for Wald intervals. Default is 0.95.
#' @param digits Number of significant digits. Default is 4.
#' @param \dots Ignored.
#'
#' @return \code{object}, invisibly.
#'
#' @seealso \code{\link{print.JMVL_Liang}}, \code{\link{JMVL_Liang}}
#'
#' @examples
#' dat <- sim_ehr_data[sim_ehr_data$id <= 80, ]
#' fit <- JMVL_Liang(dat)
#' summary(fit)
#'
#' @export
summary.JMVL_Liang <- function(object, level = 0.95, digits = 4, ...) {
  width <- 60
  cat("\nJMVL_Liang Joint Model -- Summary\n")
  cat(strrep("=", width), "\n")
  cat("\nNote: Standard errors are not available for this method.\n")
  cat("Consider bootstrap for inference.\n")

  .print_section("Visiting Process (gamma)",
                 .givehr_coef_table(object$gamma_hat, se = NULL, level), digits, width)

  if (!is.null(object$sigma_sq_hat))
    cat(sprintf("\n  Frailty variance: %.4f\n", object$sigma_sq_hat))

  .print_section("Longitudinal Outcome -- Fixed Effects (beta)",
                 .givehr_coef_table(object$beta_hat, se = NULL, level), digits, width)

  cat("\n")
  invisible(object)
}


# ============================================================================
# LMM_result
# ============================================================================

#' Print a Linear Mixed Model (LMM) for Longitudinal Data Result
#'
#' Prints the fixed-effects table from a linear mixed model fit.
#'
#' @param x An \code{LMM_result} object returned by \code{\link{LMM}}.
#' @param digits Number of significant digits. Default is 4.
#' @param level Confidence level for Wald intervals. Default is 0.95.
#' @param \dots Ignored.
#'
#' @return \code{x}, invisibly.
#'
#' @seealso \code{\link{summary.LMM_result}}, \code{\link{LMM}}
#'
#' @examples
#' dat <- sim_ehr_data[sim_ehr_data$id <= 100, ]
#' fit <- LMM(dat, fixed = Y ~ Z + X + time,
#'            random = ~ 1 | id, obs_indicator = "R")
#' print(fit)
#'
#' @export
print.LMM_result <- function(x, digits = 4, level = 0.95, ...) {
  cat("\nLinear Mixed Model -- Fixed Effects\n")
  cat(strrep("=", 50), "\n")

  tbl <- .lmm_coef_table(x$coefficients, level)
  print(round(tbl, digits))
  invisible(x)
}

#' Summarise a Linear Mixed Model (LMM) for Longitudinal Data Result
#'
#' Displays the fixed-effects coefficient table.
#'
#' @param object An \code{LMM_result} object.
#' @param level Confidence level for Wald intervals. Default is 0.95.
#' @param digits Number of significant digits. Default is 4.
#' @param \dots Ignored.
#'
#' @return \code{object}, invisibly.
#'
#' @seealso \code{\link{print.LMM_result}}, \code{\link{LMM}}
#'
#' @examples
#' dat <- sim_ehr_data[sim_ehr_data$id <= 100, ]
#' fit <- LMM(dat, fixed = Y ~ Z + X + time,
#'            random = ~ 1 | id, obs_indicator = "R")
#' summary(fit)
#'
#' @export
summary.LMM_result <- function(object, level = 0.95, digits = 4, ...) {
  width <- 60
  cat("\nLinear Mixed Model -- Summary\n")
  cat(strrep("=", width), "\n")

  .print_section("Fixed Effects",
                 .lmm_coef_table(object$coefficients, level), digits, width)

  cat("\n")
  invisible(object)
}


# ============================================================================
# SummaryStat_result
# ============================================================================

#' Print a Summary-Statistic Regression for Longitudinal Data (SummaryStat) Result
#'
#' Prints coefficients from the summary-statistic regression.
#'
#' @param x A \code{SummaryStat_result} object returned by
#'   \code{\link{SummaryStat}}.
#' @param digits Number of significant digits. Default is 4.
#' @param \dots Ignored.
#'
#' @return \code{x}, invisibly.
#'
#' @seealso \code{\link{summary.SummaryStat_result}},
#'   \code{\link{SummaryStat}}
#'
#' @examples
#' dat <- sim_ehr_data[sim_ehr_data$id <= 100, ]
#' fit <- SummaryStat(dat, formula = ~ Z + X)
#' print(fit)
#'
#' @export
print.SummaryStat_result <- function(x, digits = 4, ...) {
  cat("\nSummary-Statistic Regression -- Coefficients\n")
  cat(strrep("=", 50), "\n")

  print(round(x$coefficients, digits))
  invisible(x)
}

#' Summarise a Summary-Statistic Regression for Longitudinal Data (SummaryStat) Result
#'
#' Displays the regression coefficient table.
#'
#' @param object A \code{SummaryStat_result} object.
#' @param digits Number of significant digits. Default is 4.
#' @param \dots Ignored.
#'
#' @return \code{object}, invisibly.
#'
#' @seealso \code{\link{print.SummaryStat_result}},
#'   \code{\link{SummaryStat}}
#'
#' @examples
#' dat <- sim_ehr_data[sim_ehr_data$id <= 100, ]
#' fit <- SummaryStat(dat, formula = ~ Z + X)
#' summary(fit)
#'
#' @export
summary.SummaryStat_result <- function(object, digits = 4, ...) {
  width <- 60
  cat("\nSummary-Statistic Regression -- Summary\n")
  cat(strrep("=", width), "\n")

  cat("\nCoefficients:\n")
  cat(strrep("-", width), "\n")
  print(round(object$coefficients, digits))

  cat("\n")
  invisible(object)
}


# ============================================================================
# PairwiseLikelihood
# ============================================================================

#' Print a Pairwise Composite Likelihood (PCL) Generalized Linear Model (GLM)
#' Estimation (PairwiseLikelihood) Fit
#'
#' Prints coefficients with naive standard errors (from the inverse
#' Hessian) when available.
#'
#' @param x A \code{PairwiseLikelihood} object returned by
#'   \code{\link{PairwiseLikelihood}}.
#' @param digits Number of significant digits. Default is 4.
#' @param level Confidence level for Wald intervals. Default is 0.95.
#' @param \dots Ignored.
#'
#' @return \code{x}, invisibly.
#'
#' @seealso \code{\link{summary.PairwiseLikelihood}},
#'   \code{\link{PairwiseLikelihood}}
#'
#' @examples
#' dat <- sim_ehr_data[sim_ehr_data$id <= 100 &
#'                     sim_ehr_data$R == 1 & !is.na(sim_ehr_data$R), ]
#' fit <- PairwiseLikelihood(dat, formula = Y ~ Z + X + time,
#'                           pair_sample = 200)
#' print(fit)
#'
#' @export
print.PairwiseLikelihood <- function(x, digits = 4, level = 0.95, ...) {
  cat("\nPairwise Likelihood -- Coefficients\n")
  cat(strrep("=", 50), "\n")

  se <- NULL
  if (!is.null(x$vcov_naive))
    se <- sqrt(pmax(diag(x$vcov_naive), 0))

  tbl <- .givehr_coef_table(x$coef, se, level = level)
  print(round(tbl, digits))

  cat(sprintf("\nConverged: %s | Log-likelihood: %.2f\n",
              x$converged, x$loglik))
  invisible(x)
}

#' Summarise a Pairwise Composite Likelihood (PCL) Generalized Linear Model (GLM)
#' Estimation (PairwiseLikelihood) Fit
#'
#' Displays model dimensions, convergence status, and the coefficient
#' table with standard errors when available.
#'
#' @param object A \code{PairwiseLikelihood} object.
#' @param level Confidence level for Wald intervals. Default is 0.95.
#' @param digits Number of significant digits. Default is 4.
#' @param \dots Ignored.
#'
#' @return \code{object}, invisibly.
#'
#' @seealso \code{\link{print.PairwiseLikelihood}},
#'   \code{\link{PairwiseLikelihood}}
#'
#' @examples
#' dat <- sim_ehr_data[sim_ehr_data$id <= 100 &
#'                     sim_ehr_data$R == 1 & !is.na(sim_ehr_data$R), ]
#' fit <- PairwiseLikelihood(dat, formula = Y ~ Z + X + time,
#'                           pair_sample = 200)
#' summary(fit)
#'
#' @export
summary.PairwiseLikelihood <- function(object, level = 0.95, digits = 4, ...) {
  width <- 60
  cat("\nPairwise Likelihood -- Summary\n")
  cat(strrep("=", width), "\n")

  cat(sprintf("\n  Observations: %d | Clusters: %d | Pairs: %d\n",
              object$n, object$G, object$k_pairs))
  cat(sprintf("  Converged: %s | Log-likelihood: %.2f\n",
              object$converged, object$loglik))

  se <- NULL
  if (!is.null(object$vcov_naive))
    se <- sqrt(pmax(diag(object$vcov_naive), 0))

  .print_section("Coefficients",
                 .givehr_coef_table(object$coef, se, level), digits, width)

  cat("\n")
  invisible(object)
}


# ============================================================================
# LI_IP
# ============================================================================

#' Print a Linear Increment (LI) Imputation with Informative Presence (IP) Adjustment Fit
#'
#' Prints outcome coefficients from the linear increment imputation with
#' informative presence adjustment. Standard errors and confidence
#' intervals are not reported; see the bootstrap note in the output.
#'
#' @param x An \code{LI_IP} object returned by \code{\link{LI_IP}}.
#' @param digits Number of significant digits. Default is 4.
#' @param \dots Ignored.
#'
#' @return \code{x}, invisibly.
#'
#' @seealso \code{\link{summary.LI_IP_result}}, \code{\link{LI_IP}}
#'
#' @examples
#' dat <- sim_ehr_data[sim_ehr_data$id <= 80,
#'                     c("id", "Z", "X", "time", "R", "Y", "C")]
#' fit <- LI_IP(data = dat, ip_method = "LMM",
#'              visit_covs = c("Z", "X"))
#' print(fit)
#'
#' @export
print.LI_IP_result <- function(x, digits = 4, ...) {
  cat("\nLinear Increment + IP -- Outcome Coefficients\n")
  cat(strrep("=", 50), "\n")

  tbl <- .givehr_coef_table(x$coefficients, se = NULL)
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

#' Summarise a Linear Increment (LI) Imputation with Informative Presence (IP) Adjustment Fit
#'
#' Displays imputation model details and outcome coefficients.
#' Standard errors and confidence intervals are not reported;
#' see the bootstrap note in the output.
#'
#' @param object An \code{LI_IP} object.
#' @param digits Number of significant digits. Default is 4.
#' @param \dots Ignored.
#'
#' @return \code{object}, invisibly.
#'
#' @seealso \code{\link{print.LI_IP_result}}, \code{\link{LI_IP}}
#'
#' @examples
#' dat <- sim_ehr_data[sim_ehr_data$id <= 80,
#'                     c("id", "Z", "X", "time", "R", "Y", "C")]
#' fit <- LI_IP(data = dat, ip_method = "LMM",
#'              visit_covs = c("Z", "X"))
#' summary(fit)
#'
#' @export
summary.LI_IP_result <- function(object, digits = 4, ...) {
  width <- 60
  cat("\nLinear Increment Imputation + IP Adjustment -- Summary\n")
  cat(strrep("=", width), "\n")

  cat(sprintf("\n  LI formula: %s\n", deparse(object$li_model$formula)))
  cat(sprintf("  Covariance: %s | Values imputed: %d\n",
              object$li_model$covariance, object$n_imputed))
  cat(sprintf("  IP method: %s\n", object$ip_method))

  cat("\n  Note: Standard errors and confidence intervals are not reported.\n")
  cat("  Valid inference requires bootstrapping the full two-stage procedure\n")
  cat("  (LI imputation + IP analysis) over resampled subjects, as SE from\n")
  cat("  the IP model alone does not account for imputation uncertainty.\n")
  cat("  We suggest users write their own bootstrap method.\n")

  .print_section("Outcome Coefficients",
                 .givehr_coef_table(object$coefficients, se = NULL), digits, width)

  cat("\n")
  invisible(object)
}


# ============================================================================
# MI_IP
# ============================================================================

#' Print a Multiple Imputation (MI) with Informative Presence (IP) Adjustment Fit
#'
#' Prints pooled outcome coefficients from the multiple imputation with
#' informative presence adjustment. Standard errors and confidence
#' intervals are not reported; see the bootstrap note in the output.
#'
#' @param x An \code{MI_IP} object returned by \code{\link{MI_IP}}.
#' @param digits Number of significant digits. Default is 4.
#' @param \dots Ignored.
#'
#' @return \code{x}, invisibly.
#'
#' @seealso \code{\link{summary.MI_IP_result}}, \code{\link{MI_IP}}
#'
#' @examples
#' dat <- sim_ehr_data[sim_ehr_data$id <= 80,
#'                     c("id", "Z", "X", "time", "R", "Y", "C")]
#' fit <- MI_IP(data = dat, ip_method = "LMM",
#'              visit_covs = c("Z", "X"), m = 2, maxit = 5)
#' print(fit)
#'
#' @export
print.MI_IP_result <- function(x, digits = 4, ...) {
  cat("\nMultiple Imputation + IP -- Outcome Coefficients\n")
  cat(strrep("=", 50), "\n")

  tbl <- .givehr_coef_table(x$coefficients, se = NULL)
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

#' Summarise a Multiple Imputation (MI) with Informative Presence (IP) Adjustment Fit
#'
#' Displays imputation settings and pooled outcome coefficients.
#' Standard errors and confidence intervals are not reported;
#' see the bootstrap note in the output.
#'
#' @param object An \code{MI_IP} object.
#' @param digits Number of significant digits. Default is 4.
#' @param \dots Ignored.
#'
#' @return \code{object}, invisibly.
#'
#' @seealso \code{\link{print.MI_IP_result}}, \code{\link{MI_IP}}
#'
#' @examples
#' dat <- sim_ehr_data[sim_ehr_data$id <= 80,
#'                     c("id", "Z", "X", "time", "R", "Y", "C")]
#' fit <- MI_IP(data = dat, ip_method = "LMM",
#'              visit_covs = c("Z", "X"), m = 2, maxit = 5)
#' summary(fit)
#'
#' @export
summary.MI_IP_result <- function(object, digits = 4, ...) {
  width <- 60
  cat("\nMultiple Imputation + IP Adjustment -- Summary\n")
  cat(strrep("=", width), "\n")

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

  .print_section("Pooled Outcome Coefficients",
                 .givehr_coef_table(object$coefficients, se = NULL), digits, width)

  cat("\n")
  invisible(object)
}

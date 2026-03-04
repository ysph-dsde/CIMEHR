# ============================================================================
# Print and summary methods for benchmark / comparison methods
# ============================================================================


# ============================================================================
# EHRJoint
# ============================================================================

#' Print an EHRJoint Fit
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
#' \dontrun{
#' dat <- sim_data_gen(n = 500, seed = 42)
#' fit <- EHRJoint(dat$long_data)
#' fit
#' }
#'
#' @export
print.EHRJoint <- function(x, digits = 4, level = 0.95, ...) {
  cat("\nEHRJoint -- Outcome Fixed Effects (beta)\n")
  cat(strrep("=", 50), "\n")

  tbl <- .givehr_coef_table(x$beta_hat, se = NULL, level = level)
  print(round(tbl, digits))

  cat("\nUse summary() for all process estimates.\n")
  invisible(x)
}

#' Summarise an EHRJoint Fit
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
#' \dontrun{
#' dat <- sim_data_gen(n = 500, seed = 42)
#' fit <- EHRJoint(dat$long_data)
#' summary(fit)
#' }
#'
#' @export
summary.EHRJoint <- function(object, level = 0.95, digits = 4, ...) {
  width <- 60
  cat("\nEHRJoint Model -- Summary\n")
  cat(strrep("=", width), "\n")

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

#' Print an IIRR_weighting Fit
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
#' \dontrun{
#' dat <- sim_data_gen(n = 500, seed = 42)
#' fit <- IIRR_weighting(dat$long_data, formula = ~ Z + X)
#' fit
#' }
#'
#' @export
print.IIRR_weighting <- function(x, digits = 4, level = 0.95, ...) {
  cat("\nIIRR Weighting -- Outcome Fixed Effects (beta)\n")
  cat(strrep("=", 50), "\n")

  tbl <- .givehr_coef_table(x$beta_hat, se = NULL, level = level)
  print(round(tbl, digits))

  cat("\nUse summary() for all process estimates.\n")
  invisible(x)
}

#' Summarise an IIRR_weighting Fit
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
#' \dontrun{
#' dat <- sim_data_gen(n = 500, seed = 42)
#' fit <- IIRR_weighting(dat$long_data, formula = ~ Z + X)
#' summary(fit)
#' }
#'
#' @export
summary.IIRR_weighting <- function(object, level = 0.95, digits = 4, ...) {
  width <- 60
  cat("\nIIRR Weighting Estimator -- Summary\n")
  cat(strrep("=", width), "\n")

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

#' Print an IIRR_balancing Fit
#'
#' Prints outcome estimates from both balancing and MLE weights.
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
#' \dontrun{
#' dat <- sim_data_gen(n = 500, seed = 42)
#' fit <- IIRR_balancing(dat$long_data, formula = ~ Z + X)
#' fit
#' }
#'
#' @export
print.IIRR_balancing <- function(x, digits = 4, level = 0.95, ...) {
  cat("\nIIRR Balancing -- Outcome Fixed Effects (beta)\n")
  cat(strrep("=", 50), "\n")

  cat("\nBalancing weights:\n")
  print(round(.givehr_coef_table(x$beta_bal, se = NULL, level = level), digits))

  cat("\nMLE weights:\n")
  print(round(.givehr_coef_table(x$beta_mle, se = NULL, level = level), digits))

  cat("\nUse summary() for all process estimates.\n")
  invisible(x)
}

#' Summarise an IIRR_balancing Fit
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
#' \dontrun{
#' dat <- sim_data_gen(n = 500, seed = 42)
#' fit <- IIRR_balancing(dat$long_data, formula = ~ Z + X)
#' summary(fit)
#' }
#'
#' @export
summary.IIRR_balancing <- function(object, level = 0.95, digits = 4, ...) {
  width <- 60
  cat("\nIIRR Balancing-Weights Estimator -- Summary\n")
  cat(strrep("=", width), "\n")

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
# JMVL_IP
# ============================================================================

#' Print a JMVL_IP Fit
#'
#' Prints outcome fixed-effect estimates from the Liang--Lu--Ying joint
#' estimating equations approach.
#'
#' @param x A \code{JMVL_IP} object.
#' @param digits Number of significant digits. Default is 4.
#' @param level Confidence level for Wald intervals. Default is 0.95.
#' @param \dots Ignored.
#'
#' @return \code{x}, invisibly.
#'
#' @seealso \code{\link{summary.JMVL_IP}}, \code{\link{JMVL_IP}}
#'
#' @examples
#' \dontrun{
#' dat <- sim_data_gen(n = 500, seed = 42)
#' fit <- JMVL_IP(dat$long_data)
#' fit
#' }
#'
#' @export
print.JMVL_IP <- function(x, digits = 4, level = 0.95, ...) {
  cat("\nJMVL_IP -- Outcome Fixed Effects (beta)\n")
  cat(strrep("=", 50), "\n")

  tbl <- .givehr_coef_table(x$beta_hat, se = NULL, level = level)
  print(round(tbl, digits))

  cat("\nUse summary() for all process estimates.\n")
  invisible(x)
}

#' Summarise a JMVL_IP Fit
#'
#' Displays visiting-process and outcome estimates.
#'
#' @param object A \code{JMVL_IP} object.
#' @param level Confidence level for Wald intervals. Default is 0.95.
#' @param digits Number of significant digits. Default is 4.
#' @param \dots Ignored.
#'
#' @return \code{object}, invisibly.
#'
#' @seealso \code{\link{print.JMVL_IP}}, \code{\link{JMVL_IP}}
#'
#' @examples
#' \dontrun{
#' dat <- sim_data_gen(n = 500, seed = 42)
#' fit <- JMVL_IP(dat$long_data)
#' summary(fit)
#' }
#'
#' @export
summary.JMVL_IP <- function(object, level = 0.95, digits = 4, ...) {
  width <- 60
  cat("\nJMVL_IP Joint Model -- Summary\n")
  cat(strrep("=", width), "\n")

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

#' Print an LMM Result
#'
#' Prints the fixed-effects table from a linear mixed model fit.
#'
#' @param x An \code{LMM_result} object returned by \code{\link{LMM}}.
#' @param digits Number of significant digits. Default is 4.
#' @param \dots Ignored.
#'
#' @return \code{x}, invisibly.
#'
#' @seealso \code{\link{summary.LMM_result}}, \code{\link{LMM}}
#'
#' @examples
#' \dontrun{
#' dat <- sim_data_gen(n = 500, seed = 42)
#' fit <- LMM(dat$long_data, fixed = Y ~ Z + X + time,
#'            random = ~ 1 | id, obs_indicator = "R")
#' fit
#' }
#'
#' @export
print.LMM_result <- function(x, digits = 4, ...) {
  cat("\nLinear Mixed Model -- Fixed Effects\n")
  cat(strrep("=", 50), "\n")

  print(round(x$coefficients, digits))
  invisible(x)
}

#' Summarise an LMM Result
#'
#' Displays the fixed-effects coefficient table.
#'
#' @param object An \code{LMM_result} object.
#' @param digits Number of significant digits. Default is 4.
#' @param \dots Ignored.
#'
#' @return \code{object}, invisibly.
#'
#' @seealso \code{\link{print.LMM_result}}, \code{\link{LMM}}
#'
#' @examples
#' \dontrun{
#' dat <- sim_data_gen(n = 500, seed = 42)
#' fit <- LMM(dat$long_data, fixed = Y ~ Z + X + time,
#'            random = ~ 1 | id, obs_indicator = "R")
#' summary(fit)
#' }
#'
#' @export
summary.LMM_result <- function(object, digits = 4, ...) {
  width <- 60
  cat("\nLinear Mixed Model -- Summary\n")
  cat(strrep("=", width), "\n")

  cat("\nFixed Effects:\n")
  cat(strrep("-", width), "\n")
  print(round(object$coefficients, digits))

  cat("\n")
  invisible(object)
}


# ============================================================================
# SummaryStat_result
# ============================================================================

#' Print a SummaryStat Result
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
#' \dontrun{
#' dat <- sim_data_gen(n = 500, seed = 42)
#' fit <- SummaryStat(dat$long_data, formula = ~ Z + X)
#' fit
#' }
#'
#' @export
print.SummaryStat_result <- function(x, digits = 4, ...) {
  cat("\nSummary-Statistic Regression -- Coefficients\n")
  cat(strrep("=", 50), "\n")

  print(round(x$coefficients, digits))
  invisible(x)
}

#' Summarise a SummaryStat Result
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
#' \dontrun{
#' dat <- sim_data_gen(n = 500, seed = 42)
#' fit <- SummaryStat(dat$long_data, formula = ~ Z + X)
#' summary(fit)
#' }
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

#' Print a PairwiseLikelihood Fit
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
#' \dontrun{
#' dat <- sim_data_gen(n = 500, seed = 42)
#' fit <- PairwiseLikelihood(Y ~ Z + X, data = dat$long_data, id = "id")
#' fit
#' }
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

#' Summarise a PairwiseLikelihood Fit
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
#' \dontrun{
#' dat <- sim_data_gen(n = 500, seed = 42)
#' fit <- PairwiseLikelihood(Y ~ Z + X, data = dat$long_data, id = "id")
#' summary(fit)
#' }
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

#' Print an LI_IP Fit
#'
#' Prints outcome coefficients from the linear increment imputation with
#' informative presence adjustment.
#'
#' @param x An \code{LI_IP} object returned by \code{\link{LI_IP}}.
#' @param digits Number of significant digits. Default is 4.
#' @param level Confidence level for Wald intervals. Default is 0.95.
#' @param \dots Ignored.
#'
#' @return \code{x}, invisibly.
#'
#' @seealso \code{\link{summary.LI_IP}}, \code{\link{LI_IP}}
#'
#' @examples
#' \dontrun{
#' dat <- sim_data_gen(n = 500, seed = 42)
#' fit <- LI_IP(data = dat$long_data, ip_method = "Liang",
#'              visit_covs = c("Z", "X"))
#' fit
#' }
#'
#' @export
print.LI_IP <- function(x, digits = 4, level = 0.95, ...) {
  cat("\nLinear Increment + IP -- Outcome Coefficients\n")
  cat(strrep("=", 50), "\n")

  tbl <- .givehr_coef_table(x$coefficients, x$se, level = level)
  print(round(tbl, digits))

  cat(sprintf("\nIP method: %s | Values imputed: %d\n",
              x$ip_method, x$n_imputed))
  invisible(x)
}

#' Summarise an LI_IP Fit
#'
#' Displays imputation model details and outcome coefficients with
#' standard errors when available.
#'
#' @param object An \code{LI_IP} object.
#' @param level Confidence level for Wald intervals. Default is 0.95.
#' @param digits Number of significant digits. Default is 4.
#' @param \dots Ignored.
#'
#' @return \code{object}, invisibly.
#'
#' @seealso \code{\link{print.LI_IP}}, \code{\link{LI_IP}}
#'
#' @examples
#' \dontrun{
#' dat <- sim_data_gen(n = 500, seed = 42)
#' fit <- LI_IP(data = dat$long_data, ip_method = "Liang",
#'              visit_covs = c("Z", "X"))
#' summary(fit)
#' }
#'
#' @export
summary.LI_IP <- function(object, level = 0.95, digits = 4, ...) {
  width <- 60
  cat("\nLinear Increment Imputation + IP Adjustment -- Summary\n")
  cat(strrep("=", width), "\n")

  cat(sprintf("\n  LI formula: %s\n", deparse(object$li_model$formula)))
  cat(sprintf("  Covariance: %s | Values imputed: %d\n",
              object$li_model$covariance, object$n_imputed))
  cat(sprintf("  IP method: %s\n", object$ip_method))

  .print_section("Outcome Coefficients",
                 .givehr_coef_table(object$coefficients, object$se, level), digits, width)

  cat("\n")
  invisible(object)
}


# ============================================================================
# MI_IP
# ============================================================================

#' Print an MI_IP Fit
#'
#' Prints pooled outcome coefficients from the multiple imputation with
#' informative presence adjustment.
#'
#' @param x An \code{MI_IP} object returned by \code{\link{MI_IP}}.
#' @param digits Number of significant digits. Default is 4.
#' @param level Confidence level for Wald intervals. Default is 0.95.
#' @param \dots Ignored.
#'
#' @return \code{x}, invisibly.
#'
#' @seealso \code{\link{summary.MI_IP}}, \code{\link{MI_IP}}
#'
#' @examples
#' \dontrun{
#' dat <- sim_data_gen(n = 500, seed = 42)
#' fit <- MI_IP(data = dat$long_data, ip_method = "Liang",
#'              visit_covs = c("Z", "X"), m = 5)
#' fit
#' }
#'
#' @export
print.MI_IP <- function(x, digits = 4, level = 0.95, ...) {
  cat("\nMultiple Imputation + IP -- Outcome Coefficients\n")
  cat(strrep("=", 50), "\n")

  tbl <- .givehr_coef_table(x$coefficients, x$se, level = level)
  print(round(tbl, digits))

  cat(sprintf("\nIP method: %s | Imputations: %d/%d\n",
              x$ip_method, x$m, x$m_requested))
  invisible(x)
}

#' Summarise an MI_IP Fit
#'
#' Displays imputation settings and pooled outcome coefficients with
#' standard errors when available.
#'
#' @param object An \code{MI_IP} object.
#' @param level Confidence level for Wald intervals. Default is 0.95.
#' @param digits Number of significant digits. Default is 4.
#' @param \dots Ignored.
#'
#' @return \code{object}, invisibly.
#'
#' @seealso \code{\link{print.MI_IP}}, \code{\link{MI_IP}}
#'
#' @examples
#' \dontrun{
#' dat <- sim_data_gen(n = 500, seed = 42)
#' fit <- MI_IP(data = dat$long_data, ip_method = "Liang",
#'              visit_covs = c("Z", "X"), m = 5)
#' summary(fit)
#' }
#'
#' @export
summary.MI_IP <- function(object, level = 0.95, digits = 4, ...) {
  width <- 60
  cat("\nMultiple Imputation + IP Adjustment -- Summary\n")
  cat(strrep("=", width), "\n")

  cat(sprintf("\n  Imputation method: %s\n", object$impute_method))
  cat(sprintf("  Imputations used: %d / %d requested\n",
              object$m, object$m_requested))
  cat(sprintf("  IP method: %s\n", object$ip_method))

  .print_section("Pooled Outcome Coefficients",
                 .givehr_coef_table(object$coefficients, object$se, level), digits, width)

  cat("\n")
  invisible(object)
}

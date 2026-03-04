#' Build a Basic Coefficient Table
#'
#' Internal helper that delegates to \code{.givehr_coef_table}.
#'
#' @param est Named numeric vector of point estimates.
#' @param se Numeric vector of standard errors, or \code{NULL}.
#' @param level Numeric scalar; confidence level. Default is \code{0.95}.
#'
#' @return A \code{data.frame} with estimate, SE, and CI columns.
#'
#' @keywords internal
#' @noRd
.coef_table_basic <- function(est, se = NULL, level = 0.95) {
  .givehr_coef_table(est, se = se, level = level)
}

#' Print Method for EHRJoint Fits
#'
#' @param x An object of class \code{"EHRJoint"}.
#' @param digits Integer; number of significant digits. Default is \code{4}.
#' @param ... Additional arguments (currently unused).
#'
#' @return Invisibly returns \code{x}.
#'
#' @export
print.EHRJoint <- function(x, digits = 4, ...) {
  cat("\nEHRJoint model\n")
  cat(strrep("=", 50), "\n\n")
  cat("Call:\n"); print(x$call)
  cat("\nVisiting process (gamma):\n"); print(round(x$gamma_hat, digits))
  cat("\nObservation process (alpha):\n"); print(round(x$alpha_hat, digits))
  cat("\nOutcome process (beta):\n"); print(round(x$beta_hat, digits))
  cat("\nUse summary() for formatted tables.\n")
  invisible(x)
}


#' Print Method for IIRR_weighting Fits
#'
#' @param x An object of class \code{"IIRR_weighting"}.
#' @param digits Integer; number of significant digits. Default is \code{4}.
#' @param ... Additional arguments (currently unused).
#'
#' @return Invisibly returns \code{x}.
#'
#' @export
print.IIRR_weighting <- function(x, digits = 4, ...) {
  cat("\nIIRR weighting estimator\n")
  cat(strrep("=", 50), "\n\n")
  cat("Call:\n"); print(x$call)
  cat("\nVisiting process (gamma):\n"); print(round(x$gamma_hat, digits))
  cat("\nOutcome model (beta):\n"); print(round(x$beta_hat, digits))
  cat("\nUse summary() for formatted tables.\n")
  invisible(x)
}


#' Print Method for IIRR_balancing Fits
#'
#' @param x An object of class \code{"IIRR_balancing"}.
#' @param digits Integer; number of significant digits. Default is \code{4}.
#' @param ... Additional arguments (currently unused).
#'
#' @return Invisibly returns \code{x}.
#'
#' @export
print.IIRR_balancing <- function(x, digits = 4, ...) {
  cat("\nIIRR balancing-weights estimator\n")
  cat(strrep("=", 50), "\n\n")
  cat("Call:\n"); print(x$call)
  cat("\nVisiting process (gamma):\n"); print(round(x$gamma_hat, digits))
  cat("\nOutcome model (beta, balancing weights):\n"); print(round(x$beta_bal, digits))
  cat("\nOutcome model (beta, MLE weights):\n"); print(round(x$beta_mle, digits))
  cat("\nUse summary() for formatted tables.\n")
  invisible(x)
}


#' Print Method for JMVL_IP Fits
#'
#' @param x An object of class \code{"JMVL_IP"}.
#' @param digits Integer; number of significant digits. Default is \code{4}.
#' @param ... Additional arguments (currently unused).
#'
#' @return Invisibly returns \code{x}.
#'
#' @export
print.JMVL_IP <- function(x, digits = 4, ...) {
  cat("\nJMVL_IP joint model\n")
  cat(strrep("=", 50), "\n\n")
  cat("Call:\n"); print(x$call)
  cat("\nVisiting process (gamma):\n"); print(round(x$gamma_hat, digits))
  cat("\nOutcome process (beta):\n"); print(round(x$beta_hat, digits))
  cat("\nUse summary() for formatted tables.\n")
  invisible(x)
}

#' Print Method for LMM Results
#'
#' @param x An object of class \code{"LMM_result"}.
#' @param digits Integer; number of significant digits. Default is \code{4}.
#' @param ... Additional arguments (currently unused).
#'
#' @return Invisibly returns \code{x}.
#'
#' @export
print.LMM_result <- function(x, digits = 4, ...) {
  cat("\nLinear mixed model\n")
  cat(strrep("=", 50), "\n\n")
  cat("Call:\n"); print(x$call)
  cat("\nFixed effects table:\n")
  print(round(x$coefficients, digits))
  invisible(x)
}


#' Print Method for SummaryStat Results
#'
#' @param x An object of class \code{"SummaryStat_result"}.
#' @param digits Integer; number of significant digits. Default is \code{4}.
#' @param ... Additional arguments (currently unused).
#'
#' @return Invisibly returns \code{x}.
#'
#' @export
print.SummaryStat_result <- function(x, digits = 4, ...) {
  cat("\nSummary-statistic regression\n")
  cat(strrep("=", 50), "\n\n")
  cat("Call:\n"); print(x$call)
  cat("\nCoefficients:\n")
  print(round(x$coefficients, digits))
  invisible(x)
}

#' Print Method for PairwiseLikelihood Fits
#'
#' @param x An object of class \code{"PairwiseLikelihood"}.
#' @param digits Integer; number of significant digits. Default is \code{4}.
#' @param ... Additional arguments (currently unused).
#'
#' @return Invisibly returns \code{x}.
#'
#' @export
print.PairwiseLikelihood <- function(x, digits = 4, ...) {
  cat("\nPairwiseLikelihood fit\n")
  cat(strrep("=", 50), "\n\n")
  cat("Call:\n"); print(x$call)
  cat("\nConverged:", x$converged, "\n")
  cat("Log-likelihood:", round(x$loglik, digits), "\n\n")
  cat("Coefficients:\n")
  print(round(x$coef, digits))
  cat("\nUse summary() for SE/CI when available.\n")
  invisible(x)
}

#' Print Method for LI_IP
#'
#' @param x An object of class \code{"LI_IP"}.
#' @param ... Additional arguments (currently unused).
#'
#' @return Invisibly returns \code{x}.
#'
#' @export
print.LI_IP <- function(x, ...) {
  cat("Linear Increment Imputation with IP Adjustment\n")
  cat("===============================================\n\n")

  cat("Linear Increment Model:\n")
  cat("  Formula:", deparse(x$li_model$formula), "\n")
  cat("  Covariance:", x$li_model$covariance, "\n")
  cat("  Values imputed:", x$n_imputed, "\n\n")

  cat("IP Method:", x$ip_method, "\n")
  if (length(x$ip_args) > 0) {
    cat("IP Arguments:\n")
    for (nm in names(x$ip_args)) {
      cat("  ", nm, ":", x$ip_args[[nm]], "\n")
    }
  }
  cat("\n")

  cat("Coefficients:\n")
  if (!is.null(x$se)) {
    print_df <- data.frame(
      Estimate = x$coefficients,
      SE = x$se,
      row.names = names(x$coefficients)
    )
  } else {
    print_df <- data.frame(
      Estimate = x$coefficients,
      row.names = names(x$coefficients)
    )
  }
  print(round(print_df, 4))
  invisible(x)
}

#' Print Method for MI_IP
#'
#' @param x An object of class \code{"MI_IP_result"}.
#' @param ... Additional arguments (currently unused).
#'
#' @return Invisibly returns \code{x}.
#'
#' @export
print.MI_IP <- function(x, ...) {
  cat("Multiple Imputation with IP Adjustment\n")
  cat("=======================================\n\n")

  cat("Imputation Model:\n")
  cat("  Method:", x$impute_method, "\n")
  cat("  Fixed effects:", paste(x$impute_model$fixed, collapse = ", "), "\n")
  cat("  Random effects: intercept")
  if (!is.null(x$impute_model$random)) {
    cat(",", paste(x$impute_model$random, collapse = ", "))
  }
  cat("\n")
  cat("  Include R as predictor:", x$impute_model$include_R, "\n\n")

  cat("IP Method:", x$ip_method, "\n")
  if (length(x$ip_args) > 0) {
    cat("IP Arguments:\n")
    for (nm in names(x$ip_args)) {
      cat("  ", nm, ":", x$ip_args[[nm]], "\n")
    }
  }
  cat("Number of imputations:", x$m, "/", x$m_requested, "\n\n")

  cat("Pooled Coefficients:\n")
  if (!is.null(x$se)) {
    print_df <- data.frame(
      Estimate = x$coefficients,
      SE = x$se,
      row.names = names(x$coefficients)
    )
  } else {
    print_df <- data.frame(
      Estimate = x$coefficients,
      row.names = names(x$coefficients)
    )
  }
  print(round(print_df, 4))
  invisible(x)
}


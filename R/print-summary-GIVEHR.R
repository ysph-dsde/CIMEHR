#' Build a Coefficient Table with Optional Standard Errors and CIs
#'
#' Internal helper that constructs a data frame of estimates, standard
#' errors, and Wald confidence intervals from numeric vectors.
#'
#' @param est Named numeric vector of point estimates.
#' @param se Numeric vector of standard errors (same length as \code{est}),
#'   or \code{NULL} if unavailable.
#' @param level Numeric scalar; confidence level for Wald intervals. Default
#'   is \code{0.95}.
#'
#' @return A \code{data.frame} with columns \code{Estimate},
#'   \code{Std.Error}, \code{CI_low}, and \code{CI_high}. When \code{se}
#'   is \code{NULL}, the last three columns contain \code{NA}.
#'
#' @keywords internal
#' @noRd
.givehr_coef_table <- function(est, se = NULL, level = 0.95) {
  est <- as.numeric(est)
  nm <- names(est)
  if (!is.null(nm)) names(est) <- nm

  if (is.null(se)) {
    out <- data.frame(
      Estimate = est,
      `Std.Error` = NA_real_,
      CI_low = NA_real_,
      CI_high = NA_real_,
      check.names = FALSE
    )
    rownames(out) <- nm
    return(out)
  }

  se <- as.numeric(se)
  if (!is.null(names(se)) && is.null(nm)) nm <- names(se)
  a <- (1 - level) / 2
  z <- stats::qnorm(1 - a)

  out <- data.frame(
    Estimate = est,
    `Std.Error` = se,
    CI_low = est - z * se,
    CI_high = est + z * se,
    check.names = FALSE
  )
  rownames(out) <- nm
  out
}

#' Print Method for GIVEHR Fits
#'
#' @param x An object of class \code{"GIVEHR"}.
#' @param digits Integer; number of significant digits to print. Default is
#'   \code{4}.
#' @param ... Additional arguments (currently unused).
#'
#' @return Invisibly returns \code{x}.
#'
#' @export
print.GIVEHR <- function(x, digits = 4, ...) {
  cat("\nGIVEHR joint model\n")
  cat(strrep("=", 50), "\n\n")

  cat("Call:\n")
  print(x$call)

  cat("\nVisiting process (gamma):\n")
  print(round(x$visiting_process$gamma_V, digits))

  cat("\nObservation process (alpha):\n")
  print(round(x$observation_process$alpha_fixed, digits))

  cat("\nOutcome fixed effects (beta):\n")
  print(round(x$longitudinal_outcome$beta_Y, digits))

  cat("\nUse summary() for estimates, standard errors (when available), and confidence intervals.\n")
  invisible(x)
}

#' Print Method for GIVEHR_GH Fits
#'
#' @param x An object of class \code{"GIVEHR_GH"}.
#' @param digits Integer; number of significant digits to print. Default is
#'   \code{4}.
#' @param ... Additional arguments (currently unused).
#'
#' @return Invisibly returns \code{x}.
#'
#' @export
print.GIVEHR_GH <- function(x, digits = 4, ...) {
  cat("\nGIVEHR joint model (Gauss-Hermite)\n")
  cat(strrep("=", 50), "\n\n")

  cat("Call:\n")
  print(x$call)

  cat("\nVisiting process (gamma):\n")
  print(round(x$visiting_process$gamma_V, digits))

  cat("\nObservation process (alpha):\n")
  print(round(x$observation_process$alpha_fixed, digits))

  cat("\nOutcome fixed effects (beta):\n")
  print(round(x$longitudinal_outcome$beta_Y, digits))

  cat("\nUse summary() for a structured view of the estimates.\n")
  invisible(x)
}

# =============================================================================
# S3 Methods
# =============================================================================

#' Print Method for GIVEHR_OU_PCL Fits
#'
#' @param x An object of class \code{"GIVEHR_OU_PCL"}.
#' @param ... Additional arguments (currently unused).
#'
#' @return Invisibly returns \code{x}.
#'
#' @export
print.GIVEHR_OU_PCL <- function(x, ...) {
  cat("\nGIVEHR Model with OU+PCL Extension\n")
  cat(strrep("=", 40), "\n\n")

  cat(sprintf("Subjects: %d | Visits: %d | Observed: %d (%.1f%%)\n",
              x$model_info$n_subjects, x$model_info$n_visits, x$model_info$n_observed,
              100 * x$model_info$n_observed / x$model_info$n_visits))

  if (!is.null(x$observation_process$phi_OU) && !is.na(x$observation_process$phi_OU)) {
    rho_med <- exp(-x$observation_process$phi_OU * x$model_info$median_gap)
    cat(sprintf("OU decay (phi): %.4f | Correlation at median gap: %.4f\n",
                x$observation_process$phi_OU, rho_med))
  }

  cat("\nOutcome fixed effects (beta):\n")
  print(round(x$longitudinal_outcome$beta_Y, 4))

  cat("\nUse summary() for full results.\n")
  invisible(x)
}


#' Extract Coefficients from a GIVEHR_OU_PCL Fit
#'
#' @param object An object of class \code{"GIVEHR_OU_PCL"}.
#' @param type Character string specifying which coefficients to extract. One
#'   of \code{"outcome"} (default), \code{"visiting"}, or
#'   \code{"observation"}.
#' @param ... Additional arguments (currently unused).
#'
#' @return Named numeric vector of coefficients.
#'
#' @export
coef.GIVEHR_OU_PCL <- function(object, type = c("outcome", "visiting", "observation"), ...) {
  type <- match.arg(type)
  switch(type,
         outcome = c(object$longitudinal_outcome$beta_Y, object$longitudinal_outcome$theta_Y),
         visiting = object$visiting_process$gamma_V,
         observation = c(object$observation_process$alpha_fixed,
                         object$observation_process$delta_latent_link)
  )
}

#' Extract Variance-Covariance Matrix from a GIVEHR_OU_PCL Fit
#'
#' @param object An object of class \code{"GIVEHR_OU_PCL"}.
#' @param ... Additional arguments (currently unused).
#'
#' @return Numeric matrix; robust sandwich variance--covariance matrix for
#'   the outcome model coefficients.
#'
#' @export
vcov.GIVEHR_OU_PCL <- function(object, ...) {
  object$longitudinal_outcome$vcov_robust
}

#' Confidence Intervals for GIVEHR_OU_PCL Outcome Coefficients
#'
#' Computes Wald confidence intervals for the outcome fixed-effect
#' coefficients using the sandwich standard errors.
#'
#' @param object An object of class \code{"GIVEHR_OU_PCL"}.
#' @param parm Character vector of parameter names, or \code{NULL} (default)
#'   for all parameters.
#' @param level Numeric scalar; confidence level. Default is \code{0.95}.
#' @param ... Additional arguments (currently unused).
#'
#' @return A numeric matrix with columns for the lower and upper confidence
#'   bounds and one row per parameter.
#'
#' @export
confint.GIVEHR_OU_PCL <- function(object, parm = NULL, level = 0.95, ...) {
  beta <- object$longitudinal_outcome$beta_Y
  se <- object$longitudinal_outcome$se_beta

  a <- (1 - level) / 2
  fac <- qnorm(1 - a)

  ci <- cbind(beta - fac * se, beta + fac * se)
  colnames(ci) <- paste0(format(100 * c(a, 1 - a), trim = TRUE), "%")
  rownames(ci) <- names(beta)

  if (!is.null(parm)) ci <- ci[parm, , drop = FALSE]
  ci
}

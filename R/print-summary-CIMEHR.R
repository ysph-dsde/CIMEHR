# ============================================================================
# Print and summary methods for the CIMEHR family
#   - CIMEHR
#   - CIMEHR_timevarying_integral
#   - CIMEHR_timevarying_ou
# ============================================================================


# ============================================================================
# Internal helpers
# ============================================================================

#' Build a Coefficient Table
#'
#' @param est Named numeric vector of point estimates.
#' @param se Numeric vector of standard errors, or \code{NULL}.
#' @param level Confidence level for Wald intervals.
#'
#' @return A \code{data.frame} with column \code{Estimate} and, when
#'   \code{se} is supplied, columns \code{Std.Error}, \code{CI.low}, and
#'   \code{CI.high}.
#'
#' @keywords internal
#' @noRd
.CIMEHR_coef_table <- function(est, se = NULL, level = 0.95) {
  nm  <- names(est)
  est <- as.numeric(est)

  if (is.null(se) || all(is.na(se))) {
    out <- data.frame(Estimate = est, check.names = FALSE)
    rownames(out) <- nm
    return(out)
  }

  se <- as.numeric(se)
  z  <- stats::qnorm(1 - (1 - level) / 2)

  out <- data.frame(
    Estimate  = est,
    Std.Error = se,
    CI.low    = est - z * se,
    CI.high   = est + z * se,
    check.names = FALSE
  )
  rownames(out) <- nm
  out
}

#' Print a Section Header and Coefficient Table
#'
#' @param title Section title.
#' @param tbl A \code{data.frame} from \code{.CIMEHR_coef_table}.
#' @param digits Number of digits to round.
#' @param width Width of the separator line.
#'
#' @keywords internal
#' @noRd
.print_section <- function(title, tbl, digits = 4, width = 60) {
  cat("\n", title, "\n", strrep("-", width), "\n", sep = "")
  print(round(tbl, digits))
}

#' Validate the \code{specify} argument
#'
#' Resolves a user-supplied \code{specify} value against the layers a
#' given object actually carries. Returns a character vector of the
#' selected layer keys.
#'
#' @param specify Character; one of \code{"all"}, \code{"visiting"},
#'   \code{"observation"}, or \code{"outcome"}, or \code{NULL} (treated
#'   as \code{"all"}).
#' @param available Character vector of layer keys the object supports.
#'
#' @return Character vector of selected layer keys.
#'
#' @keywords internal
#' @noRd
.resolve_specify <- function(specify, available) {
  if (is.null(specify)) specify <- "all"
  specify <- match.arg(specify, c("all", "visiting", "observation", "outcome"))
  if (specify == "all") return(available)
  if (!specify %in% available) {
    stop(sprintf(
      "specify = \"%s\" is not available for this object. Available: %s.",
      specify, paste(sQuote(available), collapse = ", ")
    ), call. = FALSE)
  }
  specify
}


# ============================================================================
# CIMEHR
# ============================================================================

#' Print a Clinical Informative Missingness for Electronic Health Records
#' (CIMEHR) Fit
#'
#' Prints outcome fixed-effect estimates with sandwich standard errors and
#' Wald confidence intervals.
#'
#' @param x A \code{CIMEHR} object returned by \code{\link{CIMEHR}}.
#' @param digits Number of significant digits. Default is 4.
#' @param level Confidence level for Wald intervals. Default is 0.95.
#' @param \dots Ignored.
#'
#' @return \code{x}, invisibly.
#'
#' @seealso \code{\link{summary.CIMEHR}}, \code{\link{CIMEHR}}
#'
#' @examples
#' dat <- sim_ehr_data[sim_ehr_data$id <= 30, ]
#' fit <- CIMEHR(data = dat)
#' print(fit)
#'
#' @export
print.CIMEHR <- function(x, digits = 4, level = 0.95, ...) {
  cat("\nCIMEHR -- Outcome Fixed Effect (beta)\n")
  cat(strrep("=", 50), "\n")

  beta <- x$longitudinal_outcome$beta_Y
  se   <- x$longitudinal_outcome$se_beta_naive
  tbl  <- .CIMEHR_coef_table(beta, se, level = level)
  print(round(tbl, digits))

  cat("\nUse summary() for all process estimates, or\n")
  cat("summary(x, specify = \"visiting\" / \"observation\" / \"outcome\") for one layer.\n")
  invisible(x)
}

#' Summarise a Clinical Informative Missingness for Electronic Health
#' Records (CIMEHR) Fit
#'
#' Displays estimates, sandwich standard errors, and Wald confidence
#' intervals for the requested layer(s): visiting process, observation
#' process, and/or longitudinal outcome.
#'
#' @param object A \code{CIMEHR} object returned by \code{\link{CIMEHR}}.
#' @param level Confidence level for Wald intervals. Default is 0.95.
#' @param digits Number of significant digits. Default is 4.
#' @param specify Character; which process layer to display. One of
#'   \code{"all"} (default), \code{"visiting"}, \code{"observation"}, or
#'   \code{"outcome"}.
#' @param \dots Ignored.
#'
#' @return \code{object}, invisibly.
#'
#' @seealso \code{\link{print.CIMEHR}}, \code{\link{CIMEHR}}
#'
#' @examples
#' dat <- sim_ehr_data[sim_ehr_data$id <= 30, ]
#' fit <- CIMEHR(data = dat)
#' summary(fit)
#' summary(fit, specify = "visiting")
#' summary(fit, specify = "observation")
#' summary(fit, specify = "outcome")
#'
#' @export
summary.CIMEHR <- function(object, level = 0.95, digits = 4,
                           specify = c("all", "visiting", "observation", "outcome"),
                           ...) {
  specify <- match.arg(specify)
  layers  <- .resolve_specify(specify,
                              c("visiting", "observation", "outcome"))
  width   <- 60

  cat("\nCIMEHR Joint Model -- Summary\n")
  cat(strrep("=", width), "\n")
  cat("\nModel note: visiting process uses Andersen-Gill partial likelihood with log-normal frailty where estimated.\n")
  cat("Observation process uses a probit link; the OU variant adds OU serial correlation via PCL.\n")
  if (specify != "all")
    cat(sprintf("Layer requested: %s\n", specify))

  # --- Visiting process ---
  if ("visiting" %in% layers) {
    vp <- object$visiting_process
    .print_section(
      "Visiting Process -- Fixed Effect (gamma)",
      .CIMEHR_coef_table(vp$gamma_V, vp$se_gamma, level),
      digits, width
    )
    cat(sprintf("\n  Frailty variance (sigma^2_zeta): %.4f\n", vp$sigma2_zeta))
  }

  # --- Observation process ---
  if ("observation" %in% layers) {
    op <- object$observation_process
    .print_section(
      "Observation Process -- Fixed Effect (alpha)",
      .CIMEHR_coef_table(op$alpha_fixed, op$se_alpha, level),
      digits, width
    )

    if (!is.null(op$delta_latent_random_link)) {
      .print_section(
        "Observation Process -- Random Link (delta)",
        .CIMEHR_coef_table(op$delta_latent_random_link, op$se_delta, level),
        digits, width
      )
      if (!is.null(op$Sigma_q)) {
        cat("\n  Random-effect SD (Sigma_q):\n")
        print(round(op$Sigma_q, digits))
      }
    }
  }

  # --- Longitudinal outcome ---
  if ("outcome" %in% layers) {
    lo <- object$longitudinal_outcome
    beta <- lo$beta_Y
    se_beta <- lo$se_beta_naive

    .print_section(
      "Longitudinal Outcome -- Fixed Effect (beta)",
      .CIMEHR_coef_table(beta, se_beta, level),
      digits, width
    )

    if (!is.null(lo$theta_Y) && length(lo$theta_Y) > 0) {
      se_theta <- NULL
      if (!is.null(lo$vcov_joint)) {
        p_beta <- length(beta)
        idx_theta <- (p_beta + 1):ncol(lo$vcov_joint)
        se_theta <- sqrt(pmax(diag(lo$vcov_joint)[idx_theta], 0))
      }
      .print_section(
        "Longitudinal Outcome -- Random Link (theta)",
        .CIMEHR_coef_table(lo$theta_Y, se_theta, level),
        digits, width
      )
    }
  }

  cat("\n")
  invisible(object)
}


# ============================================================================
# CIMEHR_timevarying_integral
# ============================================================================

#' Print a Clinical Informative Missingness for Electronic Health Records
#' with Gauss-Hermite Quadrature (CIMEHR_timevarying_integral) Fit
#'
#' Prints outcome fixed-effect estimates from the Gauss-Hermite quadrature
#' variant. Standard errors are not available for this method.
#'
#' @param x A \code{CIMEHR_timevarying_integral} object returned by
#'   \code{\link{CIMEHR_timevarying_integral}}.
#' @param digits Number of significant digits. Default is 4.
#' @param level Confidence level for Wald intervals. Default is 0.95.
#' @param \dots Ignored.
#'
#' @return \code{x}, invisibly.
#'
#' @seealso \code{\link{summary.CIMEHR_timevarying_integral}},
#'   \code{\link{CIMEHR_timevarying_integral}}
#'
#' @examples
#' dat <- sim_ehr_data[sim_ehr_data$id <= 20, ]
#' fit <- CIMEHR_timevarying_integral(data = dat, gh_points_U = 2, gh_points_q = 2)
#' print(fit)
#'
#' @export
print.CIMEHR_timevarying_integral <- function(x, digits = 4, level = 0.95, ...) {
  cat("\nCIMEHR (Gauss-Hermite) -- Outcome Fixed Effect (beta)\n")
  cat(strrep("=", 50), "\n")

  tbl <- .CIMEHR_coef_table(x$longitudinal_outcome$beta_Y, se = NULL, level = level)
  print(round(tbl, digits))

  cat("\nNote: Standard errors are not available for CIMEHR_timevarying_integral.\n")
  cat("Consider using CIMEHR() or bootstrap for inference.\n")
  cat("Use summary() for all process estimates, or\n")
  cat("summary(x, specify = \"visiting\" / \"observation\" / \"outcome\") for one layer.\n")
  invisible(x)
}

#' Summarise a Clinical Informative Missingness for Electronic Health
#' Records with Gauss-Hermite Quadrature (CIMEHR_timevarying_integral) Fit
#'
#' Displays estimates for the requested layer(s) of the Gauss-Hermite
#' variant. Standard errors are not available for this method.
#'
#' @param object A \code{CIMEHR_timevarying_integral} object returned by
#'   \code{\link{CIMEHR_timevarying_integral}}.
#' @param level Confidence level for Wald intervals. Default is 0.95.
#' @param digits Number of significant digits. Default is 4.
#' @param specify Character; which process layer to display. One of
#'   \code{"all"} (default), \code{"visiting"}, \code{"observation"}, or
#'   \code{"outcome"}.
#' @param \dots Ignored.
#'
#' @return \code{object}, invisibly.
#'
#' @seealso \code{\link{print.CIMEHR_timevarying_integral}},
#'   \code{\link{CIMEHR_timevarying_integral}}
#'
#' @examples
#' dat <- sim_ehr_data[sim_ehr_data$id <= 20, ]
#' fit <- CIMEHR_timevarying_integral(data = dat, gh_points_U = 2, gh_points_q = 2)
#' summary(fit)
#' summary(fit, specify = "outcome")
#'
#' @export
summary.CIMEHR_timevarying_integral <- function(object, level = 0.95, digits = 4,
                                                specify = c("all", "visiting", "observation", "outcome"),
                                                ...) {
  specify <- match.arg(specify)
  layers  <- .resolve_specify(specify,
                              c("visiting", "observation", "outcome"))
  width   <- 60

  cat("\nCIMEHR Joint Model (Gauss-Hermite) -- Summary\n")
  cat(strrep("=", width), "\n")
  cat("\nModel note: visiting process uses Andersen-Gill partial likelihood with log-normal frailty where estimated.\n")
  cat("Observation process uses a probit link; the OU variant adds OU serial correlation via PCL.\n")
  if (specify != "all")
    cat(sprintf("Layer requested: %s\n", specify))
  cat("\nNote: Standard errors are not available for this method.\n")
  cat("Consider using CIMEHR() or bootstrap for inference.\n")

  if ("visiting" %in% layers) {
    vp <- object$visiting_process
    .print_section(
      "Visiting Process -- Fixed Effect (gamma)",
      .CIMEHR_coef_table(vp$gamma_V, se = NULL, level),
      digits, width
    )
    if (!is.null(vp$sigma2_zeta))
      cat(sprintf("\n  Frailty variance (sigma^2_zeta): %.4f\n", vp$sigma2_zeta))
  }

  if ("observation" %in% layers) {
    op <- object$observation_process
    .print_section(
      "Observation Process -- Fixed Effect (alpha)",
      .CIMEHR_coef_table(op$alpha_fixed, se = NULL, level),
      digits, width
    )
    if (!is.null(op$delta_latent_random_link)) {
      .print_section(
        "Observation Process -- Random Link (delta)",
        .CIMEHR_coef_table(op$delta_latent_random_link, se = NULL, level),
        digits, width
      )
    }
    if (!is.null(op$lambda))
      cat(sprintf("\n  Time-varying observation (lambda): %.4f\n", op$lambda))
  }

  if ("outcome" %in% layers) {
    lo <- object$longitudinal_outcome
    .print_section(
      "Longitudinal Outcome -- Fixed Effect (beta)",
      .CIMEHR_coef_table(lo$beta_Y, se = NULL, level),
      digits, width
    )
    if (!is.null(lo$theta_Y) && length(lo$theta_Y) > 0) {
      .print_section(
        "Longitudinal Outcome -- Random Link (theta)",
        .CIMEHR_coef_table(lo$theta_Y, se = NULL, level),
        digits, width
      )
    }
  }

  cat("\n")
  invisible(object)
}


# ============================================================================
# CIMEHR_timevarying_ou
# ============================================================================

#' Print a Clinical Informative Missingness for Electronic Health Records
#' with Ornstein-Uhlenbeck Pairwise Composite Likelihood
#' (CIMEHR_timevarying_ou) Fit
#'
#' Prints outcome fixed-effect estimates with robust standard errors from
#' the Ornstein-Uhlenbeck pairwise composite likelihood variant.
#'
#' @param x A \code{CIMEHR_timevarying_ou} object returned by
#'   \code{\link{CIMEHR_timevarying_ou}}.
#' @param digits Number of significant digits. Default is 4.
#' @param level Confidence level for Wald intervals. Default is 0.95.
#' @param \dots Ignored.
#'
#' @return \code{x}, invisibly.
#'
#' @seealso \code{\link{summary.CIMEHR_timevarying_ou}},
#'   \code{\link{CIMEHR_timevarying_ou}}
#'
#' @examples
#' dat <- sim_ehr_data[sim_ehr_data$id <= 30, ]
#' fit <- CIMEHR_timevarying_ou(data = dat, pair_type = "adjacent")
#' print(fit)
#'
#' @export
print.CIMEHR_timevarying_ou <- function(x, digits = 4, level = 0.95, ...) {
  cat("\nCIMEHR (OU + PCL) -- Outcome Fixed Effect (beta)\n")
  cat(strrep("=", 50), "\n")

  beta <- x$longitudinal_outcome$beta_Y
  se   <- x$longitudinal_outcome$se_beta
  tbl  <- .CIMEHR_coef_table(beta, se, level = level)
  print(round(tbl, digits))

  cat("\nUse summary() for all process estimates, or\n")
  cat("summary(x, specify = \"visiting\" / \"observation\" / \"outcome\") for one layer.\n")
  invisible(x)
}

#' Summarise a Clinical Informative Missingness for Electronic Health
#' Records with Ornstein-Uhlenbeck Pairwise Composite Likelihood
#' (CIMEHR_timevarying_ou) Fit
#'
#' Displays estimates for the requested layer(s) of the Ornstein-Uhlenbeck
#' variant, including the OU decay parameter. Robust standard errors are
#' shown where available.
#'
#' @param object A \code{CIMEHR_timevarying_ou} object returned by
#'   \code{\link{CIMEHR_timevarying_ou}}.
#' @param level Confidence level for Wald intervals. Default is 0.95.
#' @param digits Number of significant digits. Default is 4.
#' @param specify Character; which process layer to display. One of
#'   \code{"all"} (default), \code{"visiting"}, \code{"observation"}, or
#'   \code{"outcome"}.
#' @param \dots Ignored.
#'
#' @return \code{object}, invisibly.
#'
#' @seealso \code{\link{print.CIMEHR_timevarying_ou}},
#'   \code{\link{CIMEHR_timevarying_ou}}
#'
#' @examples
#' dat <- sim_ehr_data[sim_ehr_data$id <= 30, ]
#' fit <- CIMEHR_timevarying_ou(data = dat, pair_type = "adjacent")
#' summary(fit)
#' summary(fit, specify = "observation")
#'
#' @export
summary.CIMEHR_timevarying_ou <- function(object, level = 0.95, digits = 4,
                                          specify = c("all", "visiting", "observation", "outcome"),
                                          ...) {
  specify <- match.arg(specify)
  layers  <- .resolve_specify(specify,
                              c("visiting", "observation", "outcome"))
  width   <- 60

  cat("\nCIMEHR Joint Model (OU + PCL) -- Summary\n")
  cat(strrep("=", width), "\n")
  cat("\nModel note: visiting process uses Andersen-Gill partial likelihood with log-normal frailty where estimated.\n")
  cat("Observation process uses a probit link; the OU variant adds OU serial correlation via PCL.\n")
  if (specify != "all")
    cat(sprintf("Layer requested: %s\n", specify))

  mi <- object$model_info
  if (!is.null(mi)) {
    cat(sprintf("\n  Subjects: %d | Visits: %d | Observed: %d (%.1f%%)\n",
                mi$n_subjects, mi$n_visits, mi$n_observed,
                100 * mi$n_observed / mi$n_visits))
  }

  if ("visiting" %in% layers) {
    vp <- object$visiting_process
    .print_section(
      "Visiting Process -- Fixed Effect (gamma)",
      .CIMEHR_coef_table(vp$gamma_V, se = NULL, level),
      digits, width
    )
    if (!is.null(vp$sigma2_zeta))
      cat(sprintf("\n  Frailty variance (sigma^2_zeta): %.4f\n", vp$sigma2_zeta))
  }

  if ("observation" %in% layers) {
    op <- object$observation_process
    .print_section(
      "Observation Process -- Fixed Effect (alpha)",
      .CIMEHR_coef_table(op$alpha_fixed, se = NULL, level),
      digits, width
    )
    if (!is.null(op$delta_latent_random_link)) {
      .print_section(
        "Observation Process -- Random Link (delta)",
        .CIMEHR_coef_table(op$delta_latent_random_link, se = NULL, level),
        digits, width
      )
    }
    if (!is.null(op$phi_OU) && !is.na(op$phi_OU)) {
      cat(sprintf("\n  OU decay (phi): %.4f", op$phi_OU))
      if (!is.null(mi$median_gap))
        cat(sprintf(" | Correlation at median gap: %.4f",
                    exp(-op$phi_OU * mi$median_gap)))
      cat("\n")
    }
  }

  if ("outcome" %in% layers) {
    lo <- object$longitudinal_outcome
    .print_section(
      "Longitudinal Outcome -- Fixed Effect (beta)",
      .CIMEHR_coef_table(lo$beta_Y, lo$se_beta, level),
      digits, width
    )
    if (!is.null(lo$theta_Y) && length(lo$theta_Y) > 0) {
      se_theta <- NULL
      if (!is.null(lo$vcov_robust)) {
        p_beta <- length(lo$beta_Y)
        p_total <- ncol(lo$vcov_robust)
        if (p_total > p_beta) {
          idx_theta <- (p_beta + 1):p_total
          se_theta <- sqrt(pmax(diag(lo$vcov_robust)[idx_theta], 0))
        }
      }
      .print_section(
        "Longitudinal Outcome -- Random Link (theta)",
        .CIMEHR_coef_table(lo$theta_Y, se_theta, level),
        digits, width
      )
    }
  }

  cat("\n")
  invisible(object)
}


# ============================================================================
# coef methods
# ============================================================================

#' Extract Coefficients from a Clinical Informative Missingness for
#' Electronic Health Records (CIMEHR) Fit
#'
#' @param object A \code{CIMEHR}, \code{CIMEHR_timevarying_integral}, or
#'   \code{CIMEHR_timevarying_ou} object.
#' @param type Which coefficients: \code{"outcome"} (default),
#'   \code{"visiting"}, or \code{"observation"}.
#' @param \dots Ignored.
#'
#' @return Named numeric vector.
#'
#' @examples
#' dat <- sim_ehr_data[sim_ehr_data$id <= 30, ]
#' fit <- CIMEHR(data = dat)
#' coef(fit)
#' coef(fit, type = "visiting")
#'
#' @export
coef.CIMEHR <- function(object, type = c("outcome", "visiting", "observation"), ...) {
  type <- match.arg(type)
  switch(type,
         outcome     = object$longitudinal_outcome$beta_Y,
         visiting    = object$visiting_process$gamma_V,
         observation = object$observation_process$alpha_fixed
  )
}

#' @rdname coef.CIMEHR
#' @export
coef.CIMEHR_timevarying_integral <- function(object, type = c("outcome", "visiting", "observation"), ...) {
  type <- match.arg(type)
  switch(type,
         outcome     = object$longitudinal_outcome$beta_Y,
         visiting    = object$visiting_process$gamma_V,
         observation = object$observation_process$alpha_fixed
  )
}

#' @rdname coef.CIMEHR
#' @export
coef.CIMEHR_timevarying_ou <- function(object, type = c("outcome", "visiting", "observation"), ...) {
  type <- match.arg(type)
  switch(type,
         outcome     = c(object$longitudinal_outcome$beta_Y, object$longitudinal_outcome$theta_Y),
         visiting    = object$visiting_process$gamma_V,
         observation = c(object$observation_process$alpha_fixed,
                         object$observation_process$delta_latent_random_link)
  )
}


# ============================================================================
# vcov methods
# ============================================================================

#' Variance-Covariance Matrix for Clinical Informative Missingness for
#' Electronic Health Records (CIMEHR) Outcome Coefficients
#'
#' Returns the sandwich variance-covariance matrix for the outcome model.
#'
#' @param object A \code{CIMEHR} or \code{CIMEHR_timevarying_ou} object.
#' @param \dots Ignored.
#'
#' @return Numeric matrix.
#'
#' @examples
#' dat <- sim_ehr_data[sim_ehr_data$id <= 30, ]
#' fit <- CIMEHR(data = dat)
#' vcov(fit)
#'
#' @export
vcov.CIMEHR <- function(object, ...) {
  object$longitudinal_outcome$vcov_joint
}

#' @rdname vcov.CIMEHR
#' @export
vcov.CIMEHR_timevarying_ou <- function(object, ...) {
  object$longitudinal_outcome$vcov_robust
}


# ============================================================================
# confint methods
# ============================================================================

#' Confidence Intervals for Clinical Informative Missingness for
#' Electronic Health Records (CIMEHR) Outcome Coefficients
#'
#' Computes Wald confidence intervals using sandwich standard errors.
#'
#' @param object A \code{CIMEHR} or \code{CIMEHR_timevarying_ou} object.
#' @param parm Character vector of parameter names, or \code{NULL} for all.
#' @param level Confidence level. Default is 0.95.
#' @param \dots Ignored.
#'
#' @return Numeric matrix with lower and upper bounds, or \code{NULL} if
#'   standard errors are unavailable.
#'
#' @examples
#' dat <- sim_ehr_data[sim_ehr_data$id <= 30, ]
#' fit <- CIMEHR(data = dat)
#' confint(fit)
#' confint(fit, parm = "Z", level = 0.90)
#'
#' @export
confint.CIMEHR <- function(object, parm = NULL, level = 0.95, ...) {
  beta <- object$longitudinal_outcome$beta_Y
  se   <- object$longitudinal_outcome$se_beta_naive
  if (is.null(se)) return(NULL)

  z  <- stats::qnorm(1 - (1 - level) / 2)
  ci <- cbind(beta - z * se, beta + z * se)
  a  <- (1 - level) / 2
  colnames(ci) <- paste0(format(100 * c(a, 1 - a), trim = TRUE), "%")
  rownames(ci) <- names(beta)
  if (!is.null(parm)) ci <- ci[parm, , drop = FALSE]
  ci
}

#' @rdname confint.CIMEHR
#' @export
confint.CIMEHR_timevarying_ou <- function(object, parm = NULL, level = 0.95, ...) {
  beta <- object$longitudinal_outcome$beta_Y
  se   <- object$longitudinal_outcome$se_beta
  if (is.null(se)) return(NULL)

  z  <- stats::qnorm(1 - (1 - level) / 2)
  ci <- cbind(beta - z * se, beta + z * se)
  a  <- (1 - level) / 2
  colnames(ci) <- paste0(format(100 * c(a, 1 - a), trim = TRUE), "%")
  rownames(ci) <- names(beta)
  if (!is.null(parm)) ci <- ci[parm, , drop = FALSE]
  ci
}

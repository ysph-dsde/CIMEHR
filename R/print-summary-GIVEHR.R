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
.givehr_coef_table <- function(est, se = NULL, level = 0.95) {
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
#' @param tbl A \code{data.frame} from \code{.givehr_coef_table}.
#' @param digits Number of digits to round.
#' @param width Width of the separator line.
#'
#' @keywords internal
#' @noRd
.print_section <- function(title, tbl, digits = 4, width = 60) {
  cat("\n", title, "\n", strrep("-", width), "\n", sep = "")
  print(round(tbl, digits))
}


# ============================================================================
# GIVEHR
# ============================================================================

#' Print a GIVEHR Fit
#'
#' Prints outcome fixed-effect estimates with sandwich standard errors and
#' Wald confidence intervals.
#'
#' @param x A \code{GIVEHR} object returned by \code{\link{GIVEHR}}.
#' @param digits Number of significant digits. Default is 4.
#' @param level Confidence level for Wald intervals. Default is 0.95.
#' @param \dots Ignored.
#'
#' @return \code{x}, invisibly.
#'
#' @seealso \code{\link{summary.GIVEHR}}, \code{\link{GIVEHR}}
#'
#' @examples
#' \dontrun{
#' dat <- sim_data_gen(n = 500, seed = 42)
#' fit <- GIVEHR(data = dat$long_data)
#' fit            # outcome beta with SE and 95% CI
#' summary(fit)   # all three processes
#' }
#'
#' @export
print.GIVEHR <- function(x, digits = 4, level = 0.95, ...) {
  cat("\nGIVEHR -- Outcome Fixed Effects (beta)\n")
  cat(strrep("=", 50), "\n")

  beta <- x$longitudinal_outcome$beta_Y
  se   <- x$longitudinal_outcome$se_beta_naive
  tbl  <- .givehr_coef_table(beta, se, level = level)
  print(round(tbl, digits))

  cat("\nUse summary() for all three process estimates.\n")
  invisible(x)
}

#' Summarise a GIVEHR Fit
#'
#' Displays estimates, sandwich standard errors, and Wald confidence
#' intervals for all three stages: visiting process, observation process,
#' and longitudinal outcome.
#'
#' @param object A \code{GIVEHR} object returned by \code{\link{GIVEHR}}.
#' @param level Confidence level for Wald intervals. Default is 0.95.
#' @param digits Number of significant digits. Default is 4.
#' @param \dots Ignored.
#'
#' @return \code{object}, invisibly.
#'
#' @seealso \code{\link{print.GIVEHR}}, \code{\link{GIVEHR}}
#'
#' @examples
#' \dontrun{
#' dat <- sim_data_gen(n = 500, seed = 42)
#' fit <- GIVEHR(data = dat$long_data)
#' summary(fit)
#' summary(fit, level = 0.90)
#' }
#'
#' @export
summary.GIVEHR <- function(object, level = 0.95, digits = 4, ...) {
  width <- 60

  cat("\nGIVEHR Joint Model -- Summary\n")
  cat(strrep("=", width), "\n")

  # --- Visiting process ---
  vp <- object$visiting_process
  .print_section(
    "Visiting Process (gamma)",
    .givehr_coef_table(vp$gamma_V, vp$se_gamma, level),
    digits, width
  )
  cat(sprintf("\n  Frailty variance (sigma^2_zeta): %.4f\n", vp$sigma2_zeta))

  # --- Observation process ---
  op <- object$observation_process
  .print_section(
    "Observation Process -- Fixed Effects (alpha)",
    .givehr_coef_table(op$alpha_fixed, op$se_alpha, level),
    digits, width
  )

  if (!is.null(op$delta_latent_link)) {
    .print_section(
      "Observation Process -- Latent Link (delta)",
      .givehr_coef_table(op$delta_latent_link, op$se_delta, level),
      digits, width
    )
    if (!is.null(op$Sigma_q)) {
      cat("\n  Random-effect SD (Sigma_q):\n")
      print(round(op$Sigma_q, digits))
    }
  }

  # --- Longitudinal outcome ---
  lo <- object$longitudinal_outcome
  beta <- lo$beta_Y
  se_beta <- lo$se_beta_naive

  .print_section(
    "Longitudinal Outcome -- Fixed Effects (beta)",
    .givehr_coef_table(beta, se_beta, level),
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
      "Longitudinal Outcome -- Latent Link (theta)",
      .givehr_coef_table(lo$theta_Y, se_theta, level),
      digits, width
    )
  }

  cat("\n")
  invisible(object)
}


# ============================================================================
# GIVEHR_GH
# ============================================================================

#' Print a GIVEHR_GH Fit
#'
#' Prints outcome fixed-effect estimates from the Gauss--Hermite
#' quadrature variant.
#'
#' @param x A \code{GIVEHR_GH} object returned by \code{\link{GIVEHR_GH}}.
#' @param digits Number of significant digits. Default is 4.
#' @param level Confidence level for Wald intervals. Default is 0.95.
#' @param \dots Ignored.
#'
#' @return \code{x}, invisibly.
#'
#' @seealso \code{\link{summary.GIVEHR_GH}}, \code{\link{GIVEHR_GH}}
#'
#' @examples
#' \dontrun{
#' dat <- sim_data_gen(n = 500, seed = 42)
#' fit <- GIVEHR_GH(data = dat$long_data)
#' fit
#' }
#'
#' @export
print.GIVEHR_GH <- function(x, digits = 4, level = 0.95, ...) {
  cat("\nGIVEHR (Gauss-Hermite) -- Outcome Fixed Effects (beta)\n")
  cat(strrep("=", 50), "\n")

  tbl <- .givehr_coef_table(x$longitudinal_outcome$beta_Y, se = NULL, level = level)
  print(round(tbl, digits))

  cat("\nUse summary() for all process estimates.\n")
  invisible(x)
}

#' Summarise a GIVEHR_GH Fit
#'
#' Displays estimates for all three stages of the Gauss--Hermite variant.
#' Standard errors are shown where available.
#'
#' @param object A \code{GIVEHR_GH} object returned by
#'   \code{\link{GIVEHR_GH}}.
#' @param level Confidence level for Wald intervals. Default is 0.95.
#' @param digits Number of significant digits. Default is 4.
#' @param \dots Ignored.
#'
#' @return \code{object}, invisibly.
#'
#' @seealso \code{\link{print.GIVEHR_GH}}, \code{\link{GIVEHR_GH}}
#'
#' @examples
#' \dontrun{
#' dat <- sim_data_gen(n = 500, seed = 42)
#' fit <- GIVEHR_GH(data = dat$long_data)
#' summary(fit)
#' }
#'
#' @export
summary.GIVEHR_GH <- function(object, level = 0.95, digits = 4, ...) {
  width <- 60

  cat("\nGIVEHR Joint Model (Gauss-Hermite) -- Summary\n")
  cat(strrep("=", width), "\n")

  vp <- object$visiting_process
  .print_section(
    "Visiting Process (gamma)",
    .givehr_coef_table(vp$gamma_V, se = NULL, level),
    digits, width
  )
  if (!is.null(vp$sigma2_zeta))
    cat(sprintf("\n  Frailty variance (sigma^2_zeta): %.4f\n", vp$sigma2_zeta))

  op <- object$observation_process
  .print_section(
    "Observation Process -- Fixed Effects (alpha)",
    .givehr_coef_table(op$alpha_fixed, se = NULL, level),
    digits, width
  )
  if (!is.null(op$delta_latent_link)) {
    .print_section(
      "Observation Process -- Latent Link (delta)",
      .givehr_coef_table(op$delta_latent_link, se = NULL, level),
      digits, width
    )
  }
  if (!is.null(op$lambda))
    cat(sprintf("\n  Time-varying observation (lambda): %.4f\n", op$lambda))

  lo <- object$longitudinal_outcome
  .print_section(
    "Longitudinal Outcome -- Fixed Effects (beta)",
    .givehr_coef_table(lo$beta_Y, se = NULL, level),
    digits, width
  )
  if (!is.null(lo$theta_Y) && length(lo$theta_Y) > 0) {
    .print_section(
      "Longitudinal Outcome -- Latent Link (theta)",
      .givehr_coef_table(lo$theta_Y, se = NULL, level),
      digits, width
    )
  }

  cat("\n")
  invisible(object)
}


# ============================================================================
# GIVEHR_OU_PCL
# ============================================================================

#' Print a GIVEHR_OU_PCL Fit
#'
#' Prints outcome fixed-effect estimates with robust standard errors from
#' the OU pairwise composite likelihood variant.
#'
#' @param x A \code{GIVEHR_OU_PCL} object returned by
#'   \code{\link{GIVEHR_OU_PCL}}.
#' @param digits Number of significant digits. Default is 4.
#' @param level Confidence level for Wald intervals. Default is 0.95.
#' @param \dots Ignored.
#'
#' @return \code{x}, invisibly.
#'
#' @seealso \code{\link{summary.GIVEHR_OU_PCL}},
#'   \code{\link{GIVEHR_OU_PCL}}
#'
#' @examples
#' \dontrun{
#' dat <- sim_data_gen(n = 500, seed = 42)
#' fit <- GIVEHR_OU_PCL(data = dat$long_data)
#' fit
#' }
#'
#' @export
print.GIVEHR_OU_PCL <- function(x, digits = 4, level = 0.95, ...) {
  cat("\nGIVEHR (OU + PCL) -- Outcome Fixed Effects (beta)\n")
  cat(strrep("=", 50), "\n")

  beta <- x$longitudinal_outcome$beta_Y
  se   <- x$longitudinal_outcome$se_beta
  tbl  <- .givehr_coef_table(beta, se, level = level)
  print(round(tbl, digits))

  cat("\nUse summary() for all process estimates.\n")
  invisible(x)
}

#' Summarise a GIVEHR_OU_PCL Fit
#'
#' Displays estimates for all three stages of the OU + PCL variant,
#' including the OU decay parameter. Robust standard errors are shown
#' where available.
#'
#' @param object A \code{GIVEHR_OU_PCL} object returned by
#'   \code{\link{GIVEHR_OU_PCL}}.
#' @param level Confidence level for Wald intervals. Default is 0.95.
#' @param digits Number of significant digits. Default is 4.
#' @param \dots Ignored.
#'
#' @return \code{object}, invisibly.
#'
#' @seealso \code{\link{print.GIVEHR_OU_PCL}},
#'   \code{\link{GIVEHR_OU_PCL}}
#'
#' @examples
#' \dontrun{
#' dat <- sim_data_gen(n = 500, seed = 42)
#' fit <- GIVEHR_OU_PCL(data = dat$long_data)
#' summary(fit)
#' }
#'
#' @export
summary.GIVEHR_OU_PCL <- function(object, level = 0.95, digits = 4, ...) {
  width <- 60

  cat("\nGIVEHR Joint Model (OU + PCL) -- Summary\n")
  cat(strrep("=", width), "\n")

  mi <- object$model_info
  if (!is.null(mi)) {
    cat(sprintf("\n  Subjects: %d | Visits: %d | Observed: %d (%.1f%%)\n",
                mi$n_subjects, mi$n_visits, mi$n_observed,
                100 * mi$n_observed / mi$n_visits))
  }

  vp <- object$visiting_process
  .print_section(
    "Visiting Process (gamma)",
    .givehr_coef_table(vp$gamma_V, se = NULL, level),
    digits, width
  )
  if (!is.null(vp$sigma2_zeta))
    cat(sprintf("\n  Frailty variance (sigma^2_zeta): %.4f\n", vp$sigma2_zeta))

  op <- object$observation_process
  .print_section(
    "Observation Process -- Fixed Effects (alpha)",
    .givehr_coef_table(op$alpha_fixed, se = NULL, level),
    digits, width
  )
  if (!is.null(op$delta_latent_link)) {
    .print_section(
      "Observation Process -- Latent Link (delta)",
      .givehr_coef_table(op$delta_latent_link, se = NULL, level),
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

  lo <- object$longitudinal_outcome
  .print_section(
    "Longitudinal Outcome -- Fixed Effects (beta)",
    .givehr_coef_table(lo$beta_Y, lo$se_beta, level),
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
      "Longitudinal Outcome -- Latent Link (theta)",
      .givehr_coef_table(lo$theta_Y, se_theta, level),
      digits, width
    )
  }

  cat("\n")
  invisible(object)
}


# ============================================================================
# coef methods
# ============================================================================

#' Extract Coefficients from a GIVEHR Fit
#'
#' @param object A \code{GIVEHR}, \code{GIVEHR_GH}, or
#'   \code{GIVEHR_OU_PCL} object.
#' @param type Which coefficients: \code{"outcome"} (default),
#'   \code{"visiting"}, or \code{"observation"}.
#' @param \dots Ignored.
#'
#' @return Named numeric vector.
#'
#' @examples
#' \dontrun{
#' dat <- sim_data_gen(n = 500, seed = 42)
#' fit <- GIVEHR(data = dat$long_data)
#' coef(fit)
#' coef(fit, type = "visiting")
#' }
#'
#' @export
coef.GIVEHR <- function(object, type = c("outcome", "visiting", "observation"), ...) {
  type <- match.arg(type)
  switch(type,
         outcome     = object$longitudinal_outcome$beta_Y,
         visiting    = object$visiting_process$gamma_V,
         observation = object$observation_process$alpha_fixed
  )
}

#' @rdname coef.GIVEHR
#' @export
coef.GIVEHR_GH <- function(object, type = c("outcome", "visiting", "observation"), ...) {
  type <- match.arg(type)
  switch(type,
         outcome     = object$longitudinal_outcome$beta_Y,
         visiting    = object$visiting_process$gamma_V,
         observation = object$observation_process$alpha_fixed
  )
}

#' @rdname coef.GIVEHR
#' @export
coef.GIVEHR_OU_PCL <- function(object, type = c("outcome", "visiting", "observation"), ...) {
  type <- match.arg(type)
  switch(type,
         outcome     = c(object$longitudinal_outcome$beta_Y, object$longitudinal_outcome$theta_Y),
         visiting    = object$visiting_process$gamma_V,
         observation = c(object$observation_process$alpha_fixed,
                         object$observation_process$delta_latent_link)
  )
}


# ============================================================================
# vcov methods
# ============================================================================

#' Variance--Covariance Matrix for GIVEHR Outcome Coefficients
#'
#' Returns the sandwich variance--covariance matrix for the outcome model.
#'
#' @param object A \code{GIVEHR} or \code{GIVEHR_OU_PCL} object.
#' @param \dots Ignored.
#'
#' @return Numeric matrix.
#'
#' @examples
#' \dontrun{
#' dat <- sim_data_gen(n = 500, seed = 42)
#' fit <- GIVEHR(data = dat$long_data)
#' vcov(fit)
#' }
#'
#' @export
vcov.GIVEHR <- function(object, ...) {
  object$longitudinal_outcome$vcov_joint
}

#' @rdname vcov.GIVEHR
#' @export
vcov.GIVEHR_OU_PCL <- function(object, ...) {
  object$longitudinal_outcome$vcov_robust
}


# ============================================================================
# confint methods
# ============================================================================

#' Confidence Intervals for GIVEHR Outcome Coefficients
#'
#' Computes Wald confidence intervals using sandwich standard errors.
#'
#' @param object A \code{GIVEHR} or \code{GIVEHR_OU_PCL} object.
#' @param parm Character vector of parameter names, or \code{NULL} for all.
#' @param level Confidence level. Default is 0.95.
#' @param \dots Ignored.
#'
#' @return Numeric matrix with lower and upper bounds, or \code{NULL} if
#'   standard errors are unavailable.
#'
#' @examples
#' \dontrun{
#' dat <- sim_data_gen(n = 500, seed = 42)
#' fit <- GIVEHR(data = dat$long_data)
#' confint(fit)
#' confint(fit, parm = "Z", level = 0.90)
#' }
#'
#' @export
confint.GIVEHR <- function(object, parm = NULL, level = 0.95, ...) {
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

#' @rdname confint.GIVEHR
#' @export
confint.GIVEHR_OU_PCL <- function(object, parm = NULL, level = 0.95, ...) {
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

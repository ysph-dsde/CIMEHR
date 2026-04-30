
# ============================================================================
# Stage-specific summaries and coefficient extraction helpers
# ============================================================================

#' Summarise the Visiting Process of a Fitted Method
#'
#' @param object A fitted CIMEHR package model object.
#' @param ... Additional arguments passed to the model-specific summary method.
#'
#' @return The input object, invisibly.
#' @export
summary_visiting <- function(object, ...) {
  .summary_stage(object, stage = "visiting", ...)
}

#' Summarise the Observation Process of a Fitted Method
#'
#' @param object A fitted CIMEHR package model object.
#' @param ... Additional arguments passed to the model-specific summary method.
#'
#' @return The input object, invisibly.
#' @export
summary_observation <- function(object, ...) {
  .summary_stage(object, stage = "observation", ...)
}

#' Summarise the Outcome Process of a Fitted Method
#'
#' @param object A fitted CIMEHR package model object.
#' @param ... Additional arguments passed to the model-specific summary method.
#'
#' @return The input object, invisibly.
#' @export
summary_outcome <- function(object, ...) {
  .summary_stage(object, stage = "outcome", ...)
}

.summary_stage <- function(object, stage = c("visiting", "observation", "outcome"), ...) {
  stage <- match.arg(stage)
  if (!.method_has_stage(object, stage)) {
    cat(sprintf("\nThis method does not have a %s stage of estimation.\n", stage))
    return(invisible(object))
  }
  summary(object, specify = stage, ...)
  invisible(object)
}

.method_has_stage <- function(object, stage) {
  cls <- class(object)[1]
  stages <- switch(cls,
                   CIMEHR = c("visiting", "observation", "outcome"),
                   CIMEHR_timevarying_integral = c("visiting", "observation", "outcome"),
                   CIMEHR_timevarying_ou = c("visiting", "observation", "outcome"),
                   EHRJoint = c("visiting", "observation", "outcome"),
                   Joint_modeling_visiting_and_longitudinal_Liang = c("visiting", "outcome"),
                   Inverse_intensity_rate_ratio_weighting = c("visiting", "outcome"),
                   Inverse_intensity_rate_ratio_balancing = c("visiting", "outcome"),
                   Linear_mixed_model_result = c("outcome"),
                   Summary_stat_result = c("outcome"),
                   Pairwise_likelihood = c("outcome"),
                   Linear_increment_IP_result = c("outcome"),
                   Multiple_imputation_IP_result = c("outcome"),
                   character(0)
  )
  stage %in% stages
}

#' Extract One or More Coefficients from a Fitted Method
#'
#' This helper returns estimates from a requested model stage. Use
#' \code{parameter = NULL} to return all coefficients from that stage, or set
#' \code{parameter} to a coefficient name such as \code{"Z"}.
#'
#' @param object A fitted CIMEHR package model object.
#' @param stage One of \code{"outcome"}, \code{"visiting"}, or
#'   \code{"observation"}.
#' @param parameter Optional character vector of coefficient names to extract.
#' @param ... Ignored.
#'
#' @return A named numeric vector.
#' @export
extract_coefficient <- function(object,
                                stage = c("outcome", "visiting", "observation"),
                                parameter = NULL,
                                ...) {
  stage <- match.arg(stage)
  if (!.method_has_stage(object, stage)) {
    message(sprintf("This method does not have a %s stage of estimation.", stage))
    return(NA_real_)
  }
  est <- .extract_stage_estimates(object, stage = stage)$estimate
  if (length(est) == 0L) return(est)
  if (!is.null(parameter)) {
    missing <- setdiff(parameter, names(est))
    if (length(missing) > 0L) {
      warning("Requested parameter(s) not found: ", paste(missing, collapse = ", "), call. = FALSE)
    }
    est <- est[intersect(parameter, names(est))]
  }
  est
}

# Backward-compatible alias with a shorter name.
#' @rdname extract_coefficient
#' @export
coef_stage <- extract_coefficient

.extract_stage_estimates <- function(object, stage = c("visiting", "observation", "outcome"),
                                     which_estimate = c("balancing", "mle")) {
  stage <- match.arg(stage)
  which_estimate <- match.arg(which_estimate)
  cls <- class(object)[1]
  empty <- list(estimate = setNames(numeric(0), character(0)),
                se = NULL)

  out <- switch(cls,
                CIMEHR = .extract_stage_CIMEHR(object, stage),
                CIMEHR_timevarying_integral = .extract_stage_CIMEHR_integral(object, stage),
                CIMEHR_timevarying_ou = .extract_stage_CIMEHR_ou(object, stage),
                EHRJoint = .extract_stage_EHRJoint(object, stage),
                Joint_modeling_visiting_and_longitudinal_Liang = .extract_stage_JMVL(object, stage),
                Inverse_intensity_rate_ratio_weighting = .extract_stage_IIRR_weighting(object, stage),
                Inverse_intensity_rate_ratio_balancing = .extract_stage_IIRR_balancing(object, stage, which_estimate),
                Linear_mixed_model_result = .extract_stage_LMM(object, stage),
                Summary_stat_result = .extract_stage_SummaryStat(object, stage),
                Pairwise_likelihood = .extract_stage_Pairwise(object, stage),
                Linear_increment_IP_result = .extract_stage_LI_MI(object, stage),
                Multiple_imputation_IP_result = .extract_stage_LI_MI(object, stage),
                empty
  )
  if (is.null(out)) out <- empty
  out
}

.combine_named <- function(...) {
  x <- unlist(list(...), recursive = FALSE, use.names = TRUE)
  if (is.null(x)) setNames(numeric(0), character(0)) else unlist(x)
}

.extract_stage_CIMEHR <- function(object, stage) {
  if (stage == "visiting") return(list(estimate = object$visiting_process$gamma_V, se = object$visiting_process$se_gamma))
  if (stage == "observation") return(list(estimate = c(object$observation_process$alpha_fixed,
                                                       object$observation_process$delta_latent_random_link),
                                          se = c(object$observation_process$se_alpha,
                                                 object$observation_process$se_delta)))
  if (stage == "outcome") {
    est <- c(object$longitudinal_outcome$beta_Y, object$longitudinal_outcome$theta_Y)
    se <- object$longitudinal_outcome$se_beta_naive
    if (!is.null(object$longitudinal_outcome$vcov_joint)) {
      d <- sqrt(pmax(diag(object$longitudinal_outcome$vcov_joint), 0))
      if (length(d) == length(est)) { names(d) <- names(est); se <- d }
    }
    return(list(estimate = est, se = se))
  }
}

.extract_stage_CIMEHR_integral <- function(object, stage) {
  if (stage == "visiting") return(list(estimate = object$visiting_process$gamma_V, se = NULL))
  if (stage == "observation") return(list(estimate = c(object$observation_process$alpha_fixed,
                                                       object$observation_process$delta_latent_random_link,
                                                       lambda = object$observation_process$lambda),
                                          se = NULL))
  if (stage == "outcome") return(list(estimate = c(object$longitudinal_outcome$beta_Y,
                                                   object$longitudinal_outcome$theta_Y), se = NULL))
}

.extract_stage_CIMEHR_ou <- function(object, stage) {
  if (stage == "visiting") return(list(estimate = object$visiting_process$gamma_V, se = NULL))
  if (stage == "observation") return(list(estimate = c(object$observation_process$alpha_fixed,
                                                       object$observation_process$delta_latent_random_link,
                                                       phi_OU = object$observation_process$phi_OU),
                                          se = NULL))
  if (stage == "outcome") {
    est <- c(object$longitudinal_outcome$beta_Y, object$longitudinal_outcome$theta_Y)
    se <- object$longitudinal_outcome$se_beta
    if (!is.null(object$longitudinal_outcome$vcov_robust)) {
      d <- sqrt(pmax(diag(object$longitudinal_outcome$vcov_robust), 0))
      if (length(d) == length(est)) { names(d) <- names(est); se <- d }
    }
    return(list(estimate = est, se = se))
  }
}

.extract_stage_EHRJoint <- function(object, stage) {
  if (stage == "visiting") return(list(estimate = object$gamma_hat, se = NULL))
  if (stage == "observation") return(list(estimate = object$alpha_hat, se = NULL))
  if (stage == "outcome") return(list(estimate = c(object$beta_hat, object$theta_hat), se = NULL))
}

.extract_stage_JMVL <- function(object, stage) {
  if (stage == "visiting") return(list(estimate = object$gamma_hat, se = NULL))
  if (stage == "outcome") return(list(estimate = object$beta_hat, se = NULL))
  NULL
}

.extract_stage_IIRR_weighting <- function(object, stage) {
  if (stage == "visiting") return(list(estimate = object$gamma_hat, se = NULL))
  if (stage == "outcome") return(list(estimate = object$beta_hat, se = NULL))
  NULL
}

.extract_stage_IIRR_balancing <- function(object, stage, which_estimate = "balancing") {
  if (stage == "visiting") return(list(estimate = object$gamma_hat, se = NULL))
  if (stage == "outcome") {
    est <- if (which_estimate == "mle") object$beta_mle else object$beta_bal
    return(list(estimate = est, se = NULL))
  }
  NULL
}

.extract_stage_LMM <- function(object, stage) {
  if (stage != "outcome") return(NULL)
  df <- as.data.frame(object$coefficients, check.names = FALSE)
  est <- df[["Value"]]; names(est) <- rownames(df)
  se <- if ("Std.Error" %in% names(df)) df[["Std.Error"]] else NULL
  if (!is.null(se)) names(se) <- rownames(df)
  list(estimate = est, se = se)
}

.extract_stage_SummaryStat <- function(object, stage) {
  if (stage != "outcome") return(NULL)
  df <- as.data.frame(object$coefficients, check.names = FALSE)
  est_col <- intersect(c("Estimate", "Value"), colnames(df))[1]
  se_col <- intersect(c("Std. Error", "Std.Error"), colnames(df))[1]
  if (is.na(est_col)) return(list(estimate = setNames(numeric(0), character(0)), se = NULL))
  est <- df[[est_col]]; names(est) <- rownames(df)
  se <- if (!is.na(se_col)) df[[se_col]] else NULL
  if (!is.null(se)) names(se) <- rownames(df)
  list(estimate = est, se = se)
}

.extract_stage_Pairwise <- function(object, stage) {
  if (stage != "outcome") return(NULL)
  se <- NULL
  if (!is.null(object$vcov_naive)) se <- sqrt(pmax(diag(object$vcov_naive), 0))
  if (!is.null(se)) names(se) <- names(object$coef)
  list(estimate = object$coef, se = se)
}

.extract_stage_LI_MI <- function(object, stage) {
  if (stage != "outcome") return(NULL)
  se <- object$se
  if (!is.null(se) && is.null(names(se))) names(se) <- names(object$coefficients)
  list(estimate = object$coefficients, se = se)
}

.true_values_for_stage <- function(true_params, stage) {
  if (is.null(true_params)) return(NULL)
  if (stage == "visiting") return(true_params$gamma)
  if (stage == "observation") return(c(true_params$alpha, true_params$delta, phi_OU = true_params$phi))
  if (stage == "outcome") return(c(true_params$beta, beta_t = true_params$beta_t, true_params$theta))
  NULL
}

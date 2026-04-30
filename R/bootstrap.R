#' Bootstrap Confidence Intervals for Selected Model Fits
#'
#' Computes bootstrap standard errors and percentile confidence intervals
#' by resampling subjects with replacement and re-fitting the original model.
#' Suitable for all selected methods that do not provide closed-form standard
#' errors: \code{\link{CIMEHR_timevarying_ghq}}, \code{\link{EHRJoint}},
#' \code{\link{Joint_modeling_visiting_and_longitudinal_Liang}}, \code{\link{Inverse_intensity_rate_ratio_weighting}}, and
#' \code{\link{Inverse_intensity_rate_ratio_balancing}}. Can also supplement the sandwich standard
#' errors from \code{\link{CIMEHR}} and \code{\link{CIMEHR_timevarying_ou}}.
#'
#' For two-stage methods (\code{\link{Linear_increment_IP}}, \code{\link{Multiple_imputation_IP}}), the
#' bootstrap must cover the full two-stage procedure. We recommend users
#' implement a manual subject-level bootstrap loop.
#'
#' @param fit A fitted model object from one of the selected method functions.
#'   Supported classes: \code{"CIMEHR"}, \code{"CIMEHR_timevarying_ghq"},
#'   \code{"CIMEHR_timevarying_ou"}, \code{"EHRJoint"}, \code{"Joint_modeling_visiting_and_longitudinal_Liang"},
#'   \code{"Inverse_intensity_rate_ratio_weighting"}, \code{"Inverse_intensity_rate_ratio_balancing"}.
#' @param data The \code{data.frame} used to fit \code{fit}, in long
#'   (panel) format with one row per subject-visit.
#' @param id Character string giving the name of the subject identifier
#'   column. Default is \code{"id"}.
#' @param B Integer; number of bootstrap replicates. Default is \code{1000}.
#' @param level Confidence level for percentile intervals. Default is
#'   \code{0.95}.
#' @param covars Named list of model-specific arguments forwarded to the
#'   re-fitting function on each replicate. Supply only arguments that
#'   differ from the function defaults; unspecified arguments use the
#'   defaults listed in the Details section. Names must match the argument
#'   names of the underlying fitting function exactly; unrecognised names
#'   trigger a warning.
#' @param which_estimate For \code{"Inverse_intensity_rate_ratio_balancing"} fits only: which
#'   estimator to display in \code{print} and \code{summary}. One of
#'   \code{"bal"} (default, balancing weights) or \code{"mle"} (maximum likelihood
#'   estimation (MLE) weights). Both estimators are always stored in
#'   \code{boot_matrix} with column prefixes \code{"s3.bal."} and
#'   \code{"s3.mle."}.
#' @param seed Integer; random seed for reproducibility. Default is
#'   \code{NULL} (no seed set).
#' @param verbose Logical; if \code{TRUE}, a progress message is printed
#'   every 100 replicates. Default is \code{FALSE}.
#'
#' @return An object of class \code{"cimehr_bootstrap"}, a named list with:
#' \describe{
#'   \item{\code{estimates}}{Named numeric vector of all-stage point
#'     estimates from the original fit, with the same \code{"s1."} /
#'     \code{"s2."} / \code{"s3."} prefixes as \code{boot_matrix}.}
#'   \item{\code{boot_matrix}}{Numeric matrix of dimension
#'     \code{B_success} by \code{p}; each row is the full parameter
#'     vector from one bootstrap replicate, covering \emph{all} estimated
#'     stages. Columns are prefixed by stage: \code{"s1."} (visiting
#'     process), \code{"s2."} (observation process, where applicable),
#'     \code{"s3."} (longitudinal outcome). For \code{"Inverse_intensity_rate_ratio_balancing"},
#'     outcome columns are prefixed \code{"s3.bal."} and \code{"s3.mle."}.}
#'   \item{\code{se}}{Named numeric vector of bootstrap standard errors
#'     for the Stage-3 outcome coefficients of the displayed estimator
#'     (prefixes stripped). Standard errors for all stages are recoverable
#'     from \code{boot_matrix} directly.}
#'   \item{\code{ci}}{Numeric matrix with columns \code{CI.low} and
#'     \code{CI.high}; percentile intervals for the Stage-3 outcome
#'     coefficients of the displayed estimator (prefixes stripped).}
#'   \item{\code{level}}{The confidence level used.}
#'   \item{\code{B}}{Number of replicates requested.}
#'   \item{\code{B_success}}{Number of replicates that converged.}
#'   \item{\code{method}}{Character string; class of the original fit.}
#'   \item{\code{which_estimate}}{For \code{"Inverse_intensity_rate_ratio_balancing"}: which
#'     estimator is shown in print/summary (\code{"bal"} or \code{"mle"}).
#'     \code{NULL} for all other methods.}
#'   \item{\code{call}}{The matched call.}
#' }
#'
#' @details
#' Subjects are resampled with replacement. For each replicate, all rows
#' belonging to a sampled subject are included; subjects drawn multiple
#' times receive distinct integer identifiers (IDs) to avoid duplicate
#' identifiers in the resampled dataset.
#'
#' Replicates that throw an error or a convergence-related warning are
#' silently discarded. A warning is issued if more than 10 percent of
#' replicates fail.
#'
#' @section Covariate arguments by method:
#' Supply only entries in \code{covars} that differ from the defaults
#' shown below. Names must match the argument names of the underlying
#' fitting function exactly.
#'
#' \strong{CIMEHR, CIMEHR_timevarying_ghq, CIMEHR_timevarying_ou:}
#' \describe{
#'   \item{\code{covars_visit_XV}}{Default \code{c("Z", "X")}.}
#'   \item{\code{covars_outcome_fixed_XY}}{Default \code{c("Z", "X")}.}
#'   \item{\code{covars_outcome_link_ZY}}{Default \code{"Z"}.}
#'   \item{\code{covars_obs_fixed_XO}}{Default \code{c("Z", "X")}.}
#'   \item{\code{covars_obs_random_ZO}}{Default \code{"Z"}.}
#' }
#'
#' \strong{EHRJoint, Joint_modeling_visiting_and_longitudinal_Liang:}
#' \describe{
#'   \item{\code{visit_covs}}{Default \code{c("Z", "X")}.}
#'   \item{\code{long_covs}}{Default \code{c("Z", "X")}.}
#'   \item{\code{random_covs}}{Default \code{"Z"}.}
#' }
#'
#' \strong{Inverse_intensity_rate_ratio_weighting:}
#' \describe{
#'   \item{\code{visit_covs}}{Default \code{c("Z", "X")}.}
#'   \item{\code{formula}}{Default \code{NULL} (auto-constructed from
#'     \code{visit_covs} and \code{time}).}
#' }
#'
#' \strong{Inverse_intensity_rate_ratio_balancing:}
#' \describe{
#'   \item{\code{visit_covs}}{Default \code{c("Z", "X")}.}
#'   \item{\code{balance_covs}}{Default \code{NULL} (uses
#'     \code{visit_covs}).}
#'   \item{\code{formula}}{Default \code{NULL} (auto-constructed).}
#'   \item{\code{phi}}{Sensitivity parameter. Default \code{0}.}
#' }
#'
#' @section Plug-in vs. bootstrap:
#' \code{\link{CIMEHR}} and \code{\link{CIMEHR_timevarying_ou}} already compute
#' sandwich (plug-in) standard errors from per-subject score contributions
#' evaluated at the fitted parameters --- no re-fitting is needed.
#' Bootstrap is provided for those methods as an optional check.
#' For \code{\link{CIMEHR_timevarying_ghq}}, \code{\link{EHRJoint}},
#' \code{\link{Joint_modeling_visiting_and_longitudinal_Liang}}, \code{\link{Inverse_intensity_rate_ratio_weighting}}, and
#' \code{\link{Inverse_intensity_rate_ratio_balancing}}, per-subject score contributions are not
#' exposed by the fitting routines, so bootstrap is the recommended
#' approach for uncertainty quantification.
#'
#' @section Accessing standard errors for all stages:
#' \code{print} and \code{summary} display bootstrap standard errors (SEs) and confidence
#' intervals for the Stage-3 outcome coefficients only (\code{"s3.*"}
#' columns of \code{boot_matrix}).  Bootstrap SEs for all other stages are
#' stored in \code{boot_matrix} and can be extracted manually, for example:
#' \preformatted{
#' bs <- bootstrap(fit, data = dat, B = 500)
#'
#' # Stage-1 SE (visiting process)
#' apply(bs$boot_matrix[, grepl("^s1\\.", colnames(bs$boot_matrix)), drop = FALSE],
#'       2, sd)
#'
#' # Stage-2 SE (observation process, CIMEHR family)
#' apply(bs$boot_matrix[, grepl("^s2\\.", colnames(bs$boot_matrix)), drop = FALSE],
#'       2, sd)
#' }
#'
#' @seealso \code{\link{print.cimehr_bootstrap}},
#'   \code{\link{summary.cimehr_bootstrap}}
#'
#' @references
#' Bradley Efron and Robert J. Tibshirani (1994). \emph{An Introduction to
#' the Bootstrap}. Chapman & Hall/CRC.
#'
#' @examples
#' dat <- sim_ehr_data[sim_ehr_data$id <= 100, ]
#' fit <- Inverse_intensity_rate_ratio_weighting(dat, visit_covs = c("Z", "X"))
#' bs  <- bootstrap(fit, data = dat, B = 5, seed = 1)
#' print(bs)
#'
#' @importFrom stats quantile sd
#' @export
bootstrap <- function(fit,
                      data,
                      id             = "id",
                      B              = 1000,
                      level          = 0.95,
                      covars         = list(),
                      which_estimate = c("bal", "mle"),
                      seed           = NULL,
                      verbose        = FALSE) {

  call           <- match.call()
  which_estimate <- match.arg(which_estimate)

  # -----------------------------------------------------------------------
  # Validate class
  # -----------------------------------------------------------------------
  supported <- c("CIMEHR", "CIMEHR_timevarying_ghq", "CIMEHR_timevarying_ou",
                 "EHRJoint", "Joint_modeling_visiting_and_longitudinal_Liang",
                 "Inverse_intensity_rate_ratio_weighting", "Inverse_intensity_rate_ratio_balancing")
  fit_class <- class(fit)[1]

  if (!fit_class %in% supported) {
    stop(
      "bootstrap() does not support class '", fit_class, "'.\n",
      "Supported classes: ", paste(supported, collapse = ", "), ".\n",
      "For Linear_increment_IP and Multiple_imputation_IP, bootstrap the full two-stage procedure manually.",
      call. = FALSE
    )
  }

  if (!is.data.frame(data))
    stop("`data` must be a data.frame.", call. = FALSE)
  if (!id %in% names(data))
    stop("Column '", id, "' not found in `data`.", call. = FALSE)

  # Warn on unrecognised covars names so users catch typos early
  known   <- .bootstrap_known_covars(fit_class)
  unknown <- setdiff(names(covars), known)
  if (length(unknown) > 0)
    warning("Unrecognised covars name(s) for ", fit_class, ": ",
            paste(unknown, collapse = ", "), ". These will be ignored.",
            call. = FALSE)

  if (!is.null(seed)) set.seed(seed)

  # -----------------------------------------------------------------------
  # Point estimates from original fit
  # -----------------------------------------------------------------------
  orig_est <- .bootstrap_extract_coef(fit, fit_class, which_estimate)

  # -----------------------------------------------------------------------
  # Subject IDs for cluster resampling
  # -----------------------------------------------------------------------
  subjects <- unique(data[[id]])
  n_subj   <- length(subjects)

  # -----------------------------------------------------------------------
  # Bootstrap loop
  # -----------------------------------------------------------------------
  boot_list <- vector("list", B)
  n_fail    <- 0L

  for (b in seq_len(B)) {

    if (verbose && b %% 100 == 0)
      message("  Bootstrap replicate ", b, " / ", B)

    sampled_ids <- sample(subjects, size = n_subj, replace = TRUE)
    boot_data   <- .bootstrap_resample(data, id, sampled_ids)

    coef_b <- tryCatch(
      .bootstrap_refit(fit, fit_class, boot_data, id, covars),
      error   = function(e) NULL,
      warning = function(w) {
        if (grepl("converge|singular|failed", conditionMessage(w),
                  ignore.case = TRUE))
          return(NULL)
        withCallingHandlers(
          .bootstrap_refit(fit, fit_class, boot_data, id, covars),
          warning = function(w2) invokeRestart("muffleWarning")
        )
      }
    )

    if (is.null(coef_b)) n_fail <- n_fail + 1L else boot_list[[b]] <- coef_b
  }

  # -----------------------------------------------------------------------
  # Collect successful replicates
  # -----------------------------------------------------------------------
  boot_list <- Filter(Negate(is.null), boot_list)
  B_success <- length(boot_list)

  if (B_success == 0L)
    stop("All bootstrap replicates failed. Check that the model converges ",
         "on your data.", call. = FALSE)
  if (n_fail > 0.10 * B)
    warning(n_fail, " of ", B, " bootstrap replicates failed (>10%). ",
            "Results may be unreliable.", call. = FALSE)

  boot_mat <- do.call(rbind, boot_list)

  # -----------------------------------------------------------------------
  # SE and percentile CI for the displayed estimator only
  # -----------------------------------------------------------------------
  display_cols <- .bootstrap_display_cols(fit_class, which_estimate,
                                          colnames(boot_mat))
  alpha    <- (1 - level) / 2
  se_boot  <- apply(boot_mat[, display_cols, drop = FALSE], 2,
                    stats::sd, na.rm = TRUE)
  ci_raw   <- apply(boot_mat[, display_cols, drop = FALSE], 2,
                    stats::quantile, probs = c(alpha, 1 - alpha), na.rm = TRUE)
  ci_mat         <- t(ci_raw)
  colnames(ci_mat) <- c("CI.low", "CI.high")

  # Strip "bal." / "mle." prefix from display names
  names(se_boot)    <- .strip_estimator_prefix(names(se_boot))
  rownames(ci_mat)  <- .strip_estimator_prefix(rownames(ci_mat))

  structure(
    list(
      estimates      = orig_est,
      boot_matrix    = boot_mat,
      se             = se_boot,
      ci             = ci_mat,
      level          = level,
      B              = B,
      B_success      = B_success,
      method         = fit_class,
      which_estimate = if (fit_class == "Inverse_intensity_rate_ratio_balancing") which_estimate else NULL,
      call           = call
    ),
    class = "cimehr_bootstrap"
  )
}


# ============================================================================
# Internal helpers
# ============================================================================

#' Valid covars names per method class
#' @keywords internal
#' @noRd
.bootstrap_known_covars <- function(fit_class) {
  switch(fit_class,
         CIMEHR        = ,
         CIMEHR_timevarying_ghq     = ,
         CIMEHR_timevarying_ou = c("covars_visit_XV", "covars_outcome_fixed_XY",
                                   "covars_outcome_link_ZY", "covars_obs_fixed_XO",
                                   "covars_obs_random_ZO"),
         EHRJoint      = ,
         Joint_modeling_visiting_and_longitudinal_Liang    = c("visit_covs", "long_covs", "random_covs"),
         Inverse_intensity_rate_ratio_weighting = c("visit_covs", "formula"),
         Inverse_intensity_rate_ratio_balancing = c("visit_covs", "balance_covs", "formula", "phi"),
         character(0)
  )
}

#' Extract all-stage coefficients from original fit
#'
#' Returns a named vector covering every estimated parameter across all stages.
#' Parameters are prefixed by stage: \code{"s1."} (visiting process),
#' \code{"s2."} (observation process), \code{"s3."} (longitudinal outcome).
#' For \code{"Inverse_intensity_rate_ratio_balancing"} the outcome block is further prefixed
#' \code{"s3.bal."} and \code{"s3.mle."}.
#'
#' @keywords internal
#' @noRd
.bootstrap_extract_coef <- function(fit, fit_class, which_estimate) {
  switch(fit_class,

         # ------------------------------------------------------------------
         # CIMEHR family: three stages with named sub-lists
         # ------------------------------------------------------------------
         CIMEHR = ,
         CIMEHR_timevarying_ghq = ,
         CIMEHR_timevarying_ou = {
           s1 <- setNames(fit$visiting_process$gamma_V,
                          paste0("s1.", names(fit$visiting_process$gamma_V)))
           obs <- fit$observation_process
           s2_alpha <- setNames(obs$alpha_fixed,
                                paste0("s2.", names(obs$alpha_fixed)))
           s2_delta <- if (!is.null(obs$delta_latent_link))
             setNames(obs$delta_latent_link,
                      paste0("s2.delta.", names(obs$delta_latent_link)))
           else
             NULL
           # CIMEHR_timevarying_ou additionally carries phi_OU (scalar)
           s2_phi <- if (!is.null(obs$phi_OU))
             c(s2.phi_OU = obs$phi_OU)
           else
             NULL
           s3 <- setNames(fit$longitudinal_outcome$beta_Y,
                          paste0("s3.", names(fit$longitudinal_outcome$beta_Y)))
           c(s1, s2_alpha, s2_delta, s2_phi, s3)
         },

         # ------------------------------------------------------------------
         # EHRJoint: gamma (s1), alpha (s2), beta (s3)
         # ------------------------------------------------------------------
         EHRJoint = {
           s1 <- setNames(fit$gamma_hat,  paste0("s1.", names(fit$gamma_hat)))
           s2 <- setNames(fit$alpha_hat,  paste0("s2.", names(fit$alpha_hat)))
           s3 <- setNames(fit$beta_hat,   paste0("s3.", names(fit$beta_hat)))
           c(s1, s2, s3)
         },

         # ------------------------------------------------------------------
         # Joint_modeling_visiting_and_longitudinal_Liang: gamma (s1), beta (s3 -- no separate obs process)
         # ------------------------------------------------------------------
         Joint_modeling_visiting_and_longitudinal_Liang = {
           s1 <- setNames(fit$gamma_hat, paste0("s1.", names(fit$gamma_hat)))
           s3 <- setNames(fit$beta_hat,  paste0("s3.", names(fit$beta_hat)))
           c(s1, s3)
         },

         # ------------------------------------------------------------------
         # Inverse_intensity_rate_ratio_weighting: gamma (s1), beta (s3)
         # ------------------------------------------------------------------
         Inverse_intensity_rate_ratio_weighting = {
           s1 <- setNames(fit$gamma_hat, paste0("s1.", names(fit$gamma_hat)))
           s3 <- setNames(fit$beta_hat,  paste0("s3.", names(fit$beta_hat)))
           c(s1, s3)
         },

         # ------------------------------------------------------------------
         # Inverse_intensity_rate_ratio_balancing: gamma (s1), beta_bal and beta_mle (s3)
         # ------------------------------------------------------------------
         Inverse_intensity_rate_ratio_balancing = {
           s1      <- setNames(fit$gamma_hat,
                               paste0("s1.", names(fit$gamma_hat)))
           s3_bal  <- setNames(fit$beta_bal,
                               paste0("s3.bal.", names(fit$beta_bal)))
           s3_mle  <- setNames(fit$beta_mle,
                               paste0("s3.mle.", names(fit$beta_mle)))
           c(s1, s3_bal, s3_mle)
         },

         stop("Unknown class: ", fit_class)
  )
}

#' Select boot_matrix columns for the displayed estimator
#'
#' For \code{print} and \code{summary}, only the Stage-3 outcome columns
#' (\code{"s3.*"}) are shown.  For \code{"Inverse_intensity_rate_ratio_balancing"} the display is
#' further restricted to \code{"s3.bal.*"} or \code{"s3.mle.*"} according
#' to \code{which_estimate}.  All columns remain in \code{boot_matrix}.
#'
#' @keywords internal
#' @noRd
.bootstrap_display_cols <- function(fit_class, which_estimate, all_cols) {
  if (fit_class == "Inverse_intensity_rate_ratio_balancing") {
    prefix <- if (which_estimate == "bal") "^s3\\.bal\\." else "^s3\\.mle\\."
    return(grep(prefix, all_cols, value = TRUE))
  }
  grep("^s3\\.", all_cols, value = TRUE)
}

#' Strip stage and estimator prefixes for display
#'
#' Removes leading \code{"s3.bal."}, \code{"s3.mle."}, or \code{"s[123]."}
#' prefixes so that printed tables show clean coefficient names.
#'
#' @keywords internal
#' @noRd
.strip_estimator_prefix <- function(x) {
  x <- sub("^s3\\.(bal|mle)\\.", "", x)   # s3.bal.foo / s3.mle.foo -> foo
  x <- sub("^s[123]\\.", "", x)            # s1.foo / s2.foo / s3.foo -> foo
  x
}

#' Resample data by subject with replacement
#'
#' Subjects drawn multiple times receive new sequential integer IDs to
#' avoid duplicate subject identifiers.
#'
#' @keywords internal
#' @noRd
.bootstrap_resample <- function(data, id, sampled_ids) {
  chunks <- lapply(seq_along(sampled_ids), function(k) {
    chunk       <- data[data[[id]] == sampled_ids[k], , drop = FALSE]
    chunk[[id]] <- k
    chunk
  })
  do.call(rbind, chunks)
}

#' Re-fit a model on a bootstrap dataset
#'
#' The \code{cv()} helper resolves each argument: the user-supplied value
#' from \code{covars} takes priority, falling back to the function default.
#'
#' @keywords internal
#' @noRd
.bootstrap_refit <- function(fit, fit_class, boot_data, id, covars) {

  cv <- function(name, default) {
    if (!is.null(covars[[name]])) covars[[name]] else default
  }

  switch(fit_class,

         # ------------------------------------------------------------------
         # CIMEHR family
         # ------------------------------------------------------------------
         CIMEHR = {
           f <- CIMEHR(
             boot_data,
             id_col                  = id,
             covars_visit_XV         = cv("covars_visit_XV",         c("Z", "X")),
             covars_outcome_fixed_XY = cv("covars_outcome_fixed_XY", c("Z", "X")),
             covars_outcome_link_ZY  = cv("covars_outcome_link_ZY",  "Z"),
             covars_obs_fixed_XO     = cv("covars_obs_fixed_XO",     c("Z", "X")),
             covars_obs_random_ZO    = cv("covars_obs_random_ZO",    "Z")
           )
           .bootstrap_extract_coef(f, "CIMEHR", NULL)
         },

         CIMEHR_timevarying_ghq = {
           f <- CIMEHR_timevarying_ghq(
             boot_data,
             id_col                  = id,
             covars_visit_XV         = cv("covars_visit_XV",         c("Z", "X")),
             covars_outcome_fixed_XY = cv("covars_outcome_fixed_XY", c("Z", "X")),
             covars_outcome_link_ZY  = cv("covars_outcome_link_ZY",  "Z"),
             covars_obs_fixed_XO     = cv("covars_obs_fixed_XO",     c("Z", "X")),
             covars_obs_random_ZO    = cv("covars_obs_random_ZO",    "Z")
           )
           .bootstrap_extract_coef(f, "CIMEHR_timevarying_ghq", NULL)
         },

         CIMEHR_timevarying_ou = {
           f <- CIMEHR_timevarying_ou(
             boot_data,
             id_col                  = id,
             covars_visit_XV         = cv("covars_visit_XV",         c("Z", "X")),
             covars_outcome_fixed_XY = cv("covars_outcome_fixed_XY", c("Z", "X")),
             covars_outcome_link_ZY  = cv("covars_outcome_link_ZY",  "Z"),
             covars_obs_fixed_XO     = cv("covars_obs_fixed_XO",     c("Z", "X")),
             covars_obs_random_ZO    = cv("covars_obs_random_ZO",    "Z")
           )
           .bootstrap_extract_coef(f, "CIMEHR_timevarying_ou", NULL)
         },

         # ------------------------------------------------------------------
         # EHRJoint / Joint_modeling_visiting_and_longitudinal_Liang
         # ------------------------------------------------------------------
         EHRJoint = {
           f <- EHRJoint(
             boot_data,
             id          = id,
             visit_covs  = cv("visit_covs",  c("Z", "X")),
             long_covs   = cv("long_covs",   c("Z", "X")),
             random_covs = cv("random_covs", "Z")
           )
           .bootstrap_extract_coef(f, "EHRJoint", NULL)
         },

         Joint_modeling_visiting_and_longitudinal_Liang = {
           f <- Joint_modeling_visiting_and_longitudinal_Liang(
             boot_data,
             id          = id,
             visit_covs  = cv("visit_covs",  c("Z", "X")),
             long_covs   = cv("long_covs",   c("Z", "X")),
             random_covs = cv("random_covs", "Z")
           )
           .bootstrap_extract_coef(f, "Joint_modeling_visiting_and_longitudinal_Liang", NULL)
         },

         # ------------------------------------------------------------------
         # IIRR methods
         # ------------------------------------------------------------------
         Inverse_intensity_rate_ratio_weighting = {
           f <- Inverse_intensity_rate_ratio_weighting(
             boot_data,
             id         = id,
             visit_covs = cv("visit_covs", c("Z", "X")),
             formula    = cv("formula",    NULL)
           )
           .bootstrap_extract_coef(f, "Inverse_intensity_rate_ratio_weighting", NULL)
         },

         Inverse_intensity_rate_ratio_balancing = {
           f <- Inverse_intensity_rate_ratio_balancing(
             boot_data,
             id           = id,
             visit_covs   = cv("visit_covs",   c("Z", "X")),
             balance_covs = cv("balance_covs", NULL),
             formula      = cv("formula",      NULL),
             phi          = cv("phi",          0)
           )
           .bootstrap_extract_coef(f, "Inverse_intensity_rate_ratio_balancing", NULL)
         },

         stop("Unknown class in .bootstrap_refit: ", fit_class)
  )
}

#' Null-coalescing operator
#' @keywords internal
#' @noRd
`%||%` <- function(x, y) if (is.null(x)) y else x


# ============================================================================
# Print and summary methods
# ============================================================================

#' Print a Bootstrap Result
#'
#' Prints bootstrap standard errors and percentile confidence intervals
#' for outcome coefficients.
#'
#' @param x A \code{cimehr_bootstrap} object from \code{\link{bootstrap}}.
#' @param digits Number of significant digits. Default is 4.
#' @param \dots Ignored.
#'
#' @return \code{x}, invisibly.
#'
#' @seealso \code{\link{summary.cimehr_bootstrap}}, \code{\link{bootstrap}}
#'
#' @examples
#' dat <- sim_ehr_data[sim_ehr_data$id <= 100, ]
#' fit <- Inverse_intensity_rate_ratio_weighting(dat, visit_covs = c("Z", "X"))
#' bs  <- bootstrap(fit, data = dat, B = 5, seed = 1)
#' print(bs)
#'
#' @export
print.cimehr_bootstrap <- function(x, digits = 4, ...) {
  header <- paste0("Bootstrap Inference -- ", x$method)
  if (!is.null(x$which_estimate))
    header <- paste0(header, " (", x$which_estimate, " weights)")

  cat("\n", header, "\n", sep = "")
  cat(strrep("=", 50), "\n")

  # Subset point estimates to the same s3 display columns as se/ci
  display_cols  <- .bootstrap_display_cols(x$method,
                                           x$which_estimate %||% "bal",
                                           names(x$estimates))
  est_display   <- x$estimates[display_cols]
  names(est_display) <- .strip_estimator_prefix(names(est_display))

  tbl <- data.frame(
    Estimate = est_display,
    Boot.SE  = x$se,
    CI.low   = x$ci[, "CI.low"],
    CI.high  = x$ci[, "CI.high"],
    check.names = FALSE
  )
  print(round(tbl, digits))

  cat(sprintf("\nPercentile CI (%.0f%%) | Replicates: %d / %d succeeded\n",
              100 * x$level, x$B_success, x$B))
  cat("Full bootstrap matrix (all stages) available in $boot_matrix.\n")

  if (!is.null(x$which_estimate)) {
    other <- if (x$which_estimate == "bal") "mle" else "bal"
    cat("Use which_estimate = \"", other,
        "\" to display the other IIRR estimator.\n", sep = "")
  }

  invisible(x)
}

#' Summarise a Bootstrap Result
#'
#' Displays bootstrap standard errors, percentile confidence intervals,
#' and a summary of the replicate distribution for outcome coefficients.
#'
#' @param object A \code{cimehr_bootstrap} object from
#'   \code{\link{bootstrap}}.
#' @param digits Number of significant digits. Default is 4.
#' @param \dots Ignored.
#'
#' @return \code{object}, invisibly.
#'
#' @seealso \code{\link{print.cimehr_bootstrap}}, \code{\link{bootstrap}}
#'
#' @examples
#' dat <- sim_ehr_data[sim_ehr_data$id <= 100, ]
#' fit <- Inverse_intensity_rate_ratio_weighting(dat, visit_covs = c("Z", "X"))
#' bs  <- bootstrap(fit, data = dat, B = 5, seed = 1)
#' summary(bs)
#'
#' @export
summary.cimehr_bootstrap <- function(object, digits = 4, ...) {
  width  <- 60
  header <- paste0("Bootstrap Inference -- ", object$method, " -- Summary")
  if (!is.null(object$which_estimate))
    header <- paste0(header, " (", object$which_estimate, " weights)")

  cat("\n", header, "\n", sep = "")
  cat(strrep("=", width), "\n")

  cat(sprintf("\n  Replicates: %d requested, %d succeeded\n",
              object$B, object$B_success))
  cat(sprintf("  Confidence level: %.0f%% (percentile method)\n",
              100 * object$level))
  if (!is.null(object$which_estimate))
    cat("  Note: boot_matrix contains all stages; outcome columns prefixed 's3.bal.*' and 's3.mle.*'.\n")

  # Subset point estimates to the same s3 display columns as se/ci
  display_cols <- .bootstrap_display_cols(
    object$method,
    object$which_estimate %||% "bal",
    names(object$estimates)
  )
  est_display          <- object$estimates[display_cols]
  names(est_display)   <- .strip_estimator_prefix(names(est_display))

  tbl <- data.frame(
    Estimate = est_display,
    Boot.SE  = object$se,
    CI.low   = object$ci[, "CI.low"],
    CI.high  = object$ci[, "CI.high"],
    check.names = FALSE
  )
  cat("\nOutcome Coefficients\n")
  cat(strrep("-", width), "\n")
  print(round(tbl, digits))

  # Replicate distribution for displayed (s3) columns only
  disp_mat            <- object$boot_matrix[, display_cols, drop = FALSE]
  colnames(disp_mat)  <- .strip_estimator_prefix(colnames(disp_mat))

  cat("\nBootstrap Distribution\n")
  cat(strrep("-", width), "\n")
  print(round(apply(disp_mat, 2, summary), digits))

  cat("\n")
  invisible(object)
}

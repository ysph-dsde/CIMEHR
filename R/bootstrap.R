#' Bootstrap Confidence Intervals for Selected Model Fits
#'
#' Computes bootstrap standard errors and percentile confidence intervals
#' by resampling subjects with replacement and re-fitting the original model.
#' Suitable for all selected methods that do not provide closed-form standard
#' errors: \code{\link{GIVEHR_GH}}, \code{\link{EHRJoint}},
#' \code{\link{JMVL_Liang}}, \code{\link{IIRR_weighting}}, and
#' \code{\link{IIRR_balancing}}. Can also supplement the sandwich standard
#' errors from \code{\link{GIVEHR}} and \code{\link{GIVEHR_OU_PCL}}.
#'
#' For two-stage methods (\code{\link{LI_IP}}, \code{\link{MI_IP}}), the
#' bootstrap must cover the full two-stage procedure. We recommend users
#' implement a manual subject-level bootstrap loop.
#'
#' @param fit A fitted model object from one of the selected method functions.
#'   Supported classes: \code{"GIVEHR"}, \code{"GIVEHR_GH"},
#'   \code{"GIVEHR_OU_PCL"}, \code{"EHRJoint"}, \code{"JMVL_Liang"},
#'   \code{"IIRR_weighting"}, \code{"IIRR_balancing"}.
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
#' @param which_estimate For \code{"IIRR_balancing"} fits only: which
#'   estimator to display in \code{print} and \code{summary}. One of
#'   \code{"bal"} (default, balancing weights) or \code{"mle"} (maximum likelihood
#'   estimation (MLE) weights). Both estimators are always stored in
#'   \code{boot_matrix} with column prefixes \code{"bal."} and \code{"mle."}.
#' @param seed Integer; random seed for reproducibility. Default is
#'   \code{NULL} (no seed set).
#' @param verbose Logical; if \code{TRUE}, a progress message is printed
#'   every 100 replicates. Default is \code{FALSE}.
#'
#' @return An object of class \code{"cimehr_bootstrap"}, a named list with:
#' \describe{
#'   \item{\code{estimates}}{Named numeric vector of point estimates from
#'     the original fit (outcome coefficients of the displayed estimator).}
#'   \item{\code{boot_matrix}}{Numeric matrix of dimension
#'     \code{B_success} by \code{p}; each row is the outcome coefficient
#'     vector from one bootstrap replicate. For \code{"IIRR_balancing"},
#'     columns are prefixed \code{"bal."} and \code{"mle."}.}
#'   \item{\code{se}}{Named numeric vector of bootstrap standard errors
#'     for the displayed estimator.}
#'   \item{\code{ci}}{Numeric matrix with columns \code{CI.low} and
#'     \code{CI.high}; percentile intervals for the displayed estimator.}
#'   \item{\code{level}}{The confidence level used.}
#'   \item{\code{B}}{Number of replicates requested.}
#'   \item{\code{B_success}}{Number of replicates that converged.}
#'   \item{\code{method}}{Character string; class of the original fit.}
#'   \item{\code{which_estimate}}{For \code{"IIRR_balancing"}: which
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
#' \strong{GIVEHR, GIVEHR_GH, GIVEHR_OU_PCL:}
#' \describe{
#'   \item{\code{covars_visit_XV}}{Default \code{c("Z", "X")}.}
#'   \item{\code{covars_outcome_fixed_XY}}{Default \code{c("Z", "X")}.}
#'   \item{\code{covars_outcome_link_ZY}}{Default \code{"Z"}.}
#'   \item{\code{covars_obs_fixed_XO}}{Default \code{c("Z", "X")}.}
#'   \item{\code{covars_obs_random_ZO}}{Default \code{"Z"}.}
#' }
#'
#' \strong{EHRJoint, JMVL_Liang:}
#' \describe{
#'   \item{\code{visit_covs}}{Default \code{c("Z", "X")}.}
#'   \item{\code{long_covs}}{Default \code{c("Z", "X")}.}
#'   \item{\code{random_covs}}{Default \code{"Z"}.}
#' }
#'
#' \strong{IIRR_weighting:}
#' \describe{
#'   \item{\code{visit_covs}}{Default \code{c("Z", "X")}.}
#'   \item{\code{formula}}{Default \code{NULL} (auto-constructed from
#'     \code{visit_covs} and \code{time}).}
#' }
#'
#' \strong{IIRR_balancing:}
#' \describe{
#'   \item{\code{visit_covs}}{Default \code{c("Z", "X")}.}
#'   \item{\code{balance_covs}}{Default \code{NULL} (uses
#'     \code{visit_covs}).}
#'   \item{\code{formula}}{Default \code{NULL} (auto-constructed).}
#'   \item{\code{phi}}{Sensitivity parameter. Default \code{0}.}
#' }
#'
#' @section Plug-in vs. bootstrap:
#' \code{\link{GIVEHR}} and \code{\link{GIVEHR_OU_PCL}} already compute
#' sandwich (plug-in) standard errors from per-subject score contributions
#' evaluated at the fitted parameters --- no re-fitting is needed.
#' Bootstrap is provided for those methods as an optional check.
#' For \code{\link{GIVEHR_GH}}, \code{\link{EHRJoint}},
#' \code{\link{JMVL_Liang}}, \code{\link{IIRR_weighting}}, and
#' \code{\link{IIRR_balancing}}, per-subject score contributions are not
#' exposed by the fitting routines, so bootstrap is the recommended
#' approach for uncertainty quantification.
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
#' fit <- IIRR_weighting(dat, visit_covs = c("Z", "X"))
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
  supported <- c("GIVEHR", "GIVEHR_GH", "GIVEHR_OU_PCL",
                 "EHRJoint", "JMVL_Liang",
                 "IIRR_weighting", "IIRR_balancing")
  fit_class <- class(fit)[1]

  if (!fit_class %in% supported) {
    stop(
      "bootstrap() does not support class '", fit_class, "'.\n",
      "Supported classes: ", paste(supported, collapse = ", "), ".\n",
      "For LI_IP and MI_IP, bootstrap the full two-stage procedure manually.",
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
      which_estimate = if (fit_class == "IIRR_balancing") which_estimate else NULL,
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
         GIVEHR        = ,
         GIVEHR_GH     = ,
         GIVEHR_OU_PCL = c("covars_visit_XV", "covars_outcome_fixed_XY",
                           "covars_outcome_link_ZY", "covars_obs_fixed_XO",
                           "covars_obs_random_ZO"),
         EHRJoint      = ,
         JMVL_Liang    = c("visit_covs", "long_covs", "random_covs"),
         IIRR_weighting = c("visit_covs", "formula"),
         IIRR_balancing = c("visit_covs", "balance_covs", "formula", "phi"),
         character(0)
  )
}

#' Extract outcome coefficients from original fit for the displayed estimator
#' @keywords internal
#' @noRd
.bootstrap_extract_coef <- function(fit, fit_class, which_estimate) {
  switch(fit_class,
         GIVEHR         = fit$longitudinal_outcome$beta_Y,
         GIVEHR_GH      = fit$longitudinal_outcome$beta_Y,
         GIVEHR_OU_PCL  = fit$longitudinal_outcome$beta_Y,
         EHRJoint       = fit$beta_hat,
         JMVL_Liang     = fit$beta_hat,
         IIRR_weighting = fit$beta_hat,
         IIRR_balancing = if (which_estimate == "bal") fit$beta_bal else fit$beta_mle,
         stop("Unknown class: ", fit_class)
  )
}

#' Select boot_matrix columns for the displayed estimator
#' @keywords internal
#' @noRd
.bootstrap_display_cols <- function(fit_class, which_estimate, all_cols) {
  if (fit_class != "IIRR_balancing") return(all_cols)
  prefix <- if (which_estimate == "bal") "^bal\\." else "^mle\\."
  grep(prefix, all_cols, value = TRUE)
}

#' Strip "bal." / "mle." prefix from coefficient names
#' @keywords internal
#' @noRd
.strip_estimator_prefix <- function(x) sub("^(bal|mle)\\.", "", x)

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
         # GIVEHR family
         # ------------------------------------------------------------------
         GIVEHR = {
           f <- GIVEHR(
             boot_data,
             id_col                  = id,
             covars_visit_XV         = cv("covars_visit_XV",         c("Z", "X")),
             covars_outcome_fixed_XY = cv("covars_outcome_fixed_XY", c("Z", "X")),
             covars_outcome_link_ZY  = cv("covars_outcome_link_ZY",  "Z"),
             covars_obs_fixed_XO     = cv("covars_obs_fixed_XO",     c("Z", "X")),
             covars_obs_random_ZO    = cv("covars_obs_random_ZO",    "Z")
           )
           f$longitudinal_outcome$beta_Y
         },

         GIVEHR_GH = {
           f <- GIVEHR_GH(
             boot_data,
             id_col                  = id,
             covars_visit_XV         = cv("covars_visit_XV",         c("Z", "X")),
             covars_outcome_fixed_XY = cv("covars_outcome_fixed_XY", c("Z", "X")),
             covars_outcome_link_ZY  = cv("covars_outcome_link_ZY",  "Z"),
             covars_obs_fixed_XO     = cv("covars_obs_fixed_XO",     c("Z", "X")),
             covars_obs_random_ZO    = cv("covars_obs_random_ZO",    "Z")
           )
           f$longitudinal_outcome$beta_Y
         },

         GIVEHR_OU_PCL = {
           f <- GIVEHR_OU_PCL(
             boot_data,
             id_col                  = id,
             covars_visit_XV         = cv("covars_visit_XV",         c("Z", "X")),
             covars_outcome_fixed_XY = cv("covars_outcome_fixed_XY", c("Z", "X")),
             covars_outcome_link_ZY  = cv("covars_outcome_link_ZY",  "Z"),
             covars_obs_fixed_XO     = cv("covars_obs_fixed_XO",     c("Z", "X")),
             covars_obs_random_ZO    = cv("covars_obs_random_ZO",    "Z")
           )
           f$longitudinal_outcome$beta_Y
         },

         # ------------------------------------------------------------------
         # EHRJoint / JMVL_Liang
         # ------------------------------------------------------------------
         EHRJoint = {
           f <- EHRJoint(
             boot_data,
             id          = id,
             visit_covs  = cv("visit_covs",  c("Z", "X")),
             long_covs   = cv("long_covs",   c("Z", "X")),
             random_covs = cv("random_covs", "Z")
           )
           f$beta_hat
         },

         JMVL_Liang = {
           f <- JMVL_Liang(
             boot_data,
             id          = id,
             visit_covs  = cv("visit_covs",  c("Z", "X")),
             long_covs   = cv("long_covs",   c("Z", "X")),
             random_covs = cv("random_covs", "Z")
           )
           f$beta_hat
         },

         # ------------------------------------------------------------------
         # IIRR methods
         # ------------------------------------------------------------------
         IIRR_weighting = {
           f <- IIRR_weighting(
             boot_data,
             id         = id,
             visit_covs = cv("visit_covs", c("Z", "X")),
             formula    = cv("formula",    NULL)
           )
           f$beta_hat
         },

         IIRR_balancing = {
           f <- IIRR_balancing(
             boot_data,
             id           = id,
             visit_covs   = cv("visit_covs",   c("Z", "X")),
             balance_covs = cv("balance_covs", NULL),
             formula      = cv("formula",      NULL),
             phi          = cv("phi",          0)
           )
           # Store both estimators; columns distinguished by prefix
           c(setNames(f$beta_bal, paste0("bal.", names(f$beta_bal))),
             setNames(f$beta_mle, paste0("mle.", names(f$beta_mle))))
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
#' fit <- IIRR_weighting(dat, visit_covs = c("Z", "X"))
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

  tbl <- data.frame(
    Estimate = x$estimates,
    Boot.SE  = x$se,
    CI.low   = x$ci[, "CI.low"],
    CI.high  = x$ci[, "CI.high"],
    check.names = FALSE
  )
  print(round(tbl, digits))

  cat(sprintf("\nPercentile CI (%.0f%%) | Replicates: %d / %d succeeded\n",
              100 * x$level, x$B_success, x$B))

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
#' fit <- IIRR_weighting(dat, visit_covs = c("Z", "X"))
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
    cat("  Note: boot_matrix contains both 'bal.*' and 'mle.*' columns.\n")

  tbl <- data.frame(
    Estimate = object$estimates,
    Boot.SE  = object$se,
    CI.low   = object$ci[, "CI.low"],
    CI.high  = object$ci[, "CI.high"],
    check.names = FALSE
  )
  cat("\nOutcome Coefficients\n")
  cat(strrep("-", width), "\n")
  print(round(tbl, digits))

  # Replicate distribution for displayed columns only
  display_cols <- .bootstrap_display_cols(
    object$method,
    object$which_estimate %||% "bal",
    colnames(object$boot_matrix)
  )
  disp_mat            <- object$boot_matrix[, display_cols, drop = FALSE]
  colnames(disp_mat)  <- .strip_estimator_prefix(colnames(disp_mat))

  cat("\nBootstrap Distribution\n")
  cat(strrep("-", width), "\n")
  print(round(apply(disp_mat, 2, summary), digits))

  cat("\n")
  invisible(object)
}

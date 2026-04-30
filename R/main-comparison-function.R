#' Available Methods for method_comparisons()
#'
#' @return Character vector of method names that can be fitted by
#'   \code{method_comparisons()}.
#' @export
available_comparison_methods <- function() {
  c(
    "CIMEHR",
    "CIMEHR_timevarying_integral",
    "CIMEHR_timevarying_ou",
    "EHRJoint",
    "Inverse_intensity_rate_ratio_balancing",
    "Inverse_intensity_rate_ratio_weighting",
    "Joint_modeling_visiting_and_longitudinal_Liang",
    "Linear_increment",
    "Linear_mixed_model",
    "Multiple_imputation",
    "Pairwise_likelihood",
    "Summary_stat"
  )
}

#' Fit and Compare Multiple Methods
#'
#' Fits any combination of implemented CIMEHR package methods and returns
#' process-specific comparison tables. Each table contains one row per method
#' and parameter. For methods that do not estimate a requested process, the
#' table contains an \code{NA} row for that process. When true values are
#' available, bias and RMSE columns are included.
#'
#' @param data Long-format analysis data. This can also be a \code{sim_CIMEHR}
#'   object returned by \code{sim_data_gen()}, in which case \code{data$long_data}
#'   is used and true simulation parameters are detected automatically.
#' @param methods Character vector of method names to fit, or \code{"all"}.
#'   Available names are returned by \code{available_comparison_methods()}.
#'   The aliases \code{"Linear_increment"} and \code{"Multiple_imputation"}
#'   call \code{Linear_increment_IP()} and \code{Multiple_imputation_IP()},
#'   respectively.
#' @param method_args Optional named list of argument lists. Names should match
#'   entries in \code{methods}. Each sub-list is passed to the corresponding
#'   fitting function. Names may use either the display name (e.g.,
#'   \code{"Linear_increment"}) or the underlying function name (e.g.,
#'   \code{"Linear_increment_IP"}).
#' @param printSE Logical; include standard error column where available.
#' @param print95CI Logical; include Wald 95 percent confidence interval columns
#'   where standard errors are available.
#' @param true_values Optional named list of true parameters. If \code{NULL},
#'   the function tries to use \code{attr(data, "true_params")} or
#'   \code{data$params} for simulated data.
#' @param true_value_verbose Logical; include true value, bias, and RMSE
#'   columns when true values are available.
#' @param which_balancing Character; for
#'   \code{Inverse_intensity_rate_ratio_balancing}, choose outcome estimates
#'   from \code{"balancing"} or \code{"mle"} weights.
#' @param default_formula Formula used for methods that require a formula when
#'   one is not supplied in \code{method_args}; default is \code{Y ~ Z + X + time}.
#' @param default_random_formula Random-effects formula used for
#'   \code{Linear_mixed_model()} and imputation methods when one is not supplied;
#'   default is \code{~ 1 | id}.
#' @param stop_on_error Logical; if \code{TRUE}, stop when a method fails. If
#'   \code{FALSE}, store the error and continue fitting the remaining methods.
#' @param ... Additional arguments passed to all methods unless overridden by
#'   method-specific entries in \code{method_args}.
#'
#' @return An object of class \code{"CIMEHR_method_comparisons"} with fitted
#'   model objects and process-specific tables.
#' @export
method_comparisons <- function(data,
                               methods = "all",
                               method_args = list(),
                               printSE = TRUE,
                               print95CI = TRUE,
                               true_values = NULL,
                               true_value_verbose = TRUE,
                               which_balancing = c("balancing", "mle"),
                               default_formula = Y ~ Z + X + time,
                               default_random_formula = ~ 1 | id,
                               stop_on_error = FALSE,
                               ...) {
  which_balancing <- match.arg(which_balancing)
  extra_args <- list(...)

  if (length(methods) == 1L && identical(tolower(methods), "all")) {
    methods <- available_comparison_methods()
  }
  methods <- as.character(methods)

  # Allow users to pass the full object returned by sim_data_gen().
  if (inherits(data, "sim_CIMEHR") && !is.null(data$long_data)) {
    if (is.null(true_values)) true_values <- data$params
    data <- data$long_data
  }
  if (is.null(true_values)) true_values <- attr(data, "true_params")

  fits <- list()
  errors <- list()

  for (display_name in methods) {
    spec <- .resolve_comparison_method(display_name)
    if (is.null(spec)) {
      msg <- paste0(
        "Method not recognized: ", display_name,
        ". Available methods are: ",
        paste(available_comparison_methods(), collapse = ", ")
      )
      if (stop_on_error) stop(msg, call. = FALSE)
      warning(msg, call. = FALSE)
      fits[[display_name]] <- NULL
      errors[[display_name]] <- msg
      next
    }

    f <- spec$fun
    if (!exists(f, mode = "function")) {
      msg <- paste0("Method function not found: ", f,
                    " (requested as ", display_name, ")")
      if (stop_on_error) stop(msg, call. = FALSE)
      warning(msg, call. = FALSE)
      fits[[display_name]] <- NULL
      errors[[display_name]] <- msg
      next
    }

    args_m <- .build_method_args(
      display_name = display_name,
      fun_name = f,
      data = data,
      extra_args = extra_args,
      method_args = method_args,
      default_formula = default_formula,
      default_random_formula = default_random_formula
    )

    fit_m <- tryCatch(do.call(get(f, mode = "function"), args_m),
                      error = function(e) e)
    if (inherits(fit_m, "error")) {
      if (stop_on_error) stop(fit_m)
      warning("Method ", display_name, " failed: ",
              conditionMessage(fit_m), call. = FALSE)
      fits[[display_name]] <- NULL
      errors[[display_name]] <- conditionMessage(fit_m)
    } else {
      fits[[display_name]] <- fit_m
    }
  }

  stages <- c("visiting", "observation", "outcome")
  tables <- setNames(vector("list", length(stages)), stages)
  for (stage in stages) {
    rows <- list()
    for (display_name in methods) {
      fit <- fits[[display_name]]
      if (is.null(fit)) {
        rows[[length(rows) + 1L]] <- .empty_comparison_row(display_name)
        next
      }
      if (!.method_has_stage(fit, stage)) {
        rows[[length(rows) + 1L]] <- .empty_comparison_row(display_name)
        next
      }

      ex <- .extract_stage_estimates(fit, stage = stage,
                                     which_estimate = which_balancing)
      est <- ex$estimate
      se <- ex$se
      if (length(est) == 0L) {
        rows[[length(rows) + 1L]] <- .empty_comparison_row(display_name)
        next
      }
      if (is.null(names(est)) || any(!nzchar(names(est)))) {
        names(est) <- paste0("param", seq_along(est))
      }

      df <- data.frame(
        Method = display_name,
        `Parameter of interest` = names(est),
        Estimate = as.numeric(est),
        check.names = FALSE,
        stringsAsFactors = FALSE
      )

      se_vec <- rep(NA_real_, length(est)); names(se_vec) <- names(est)
      if (!is.null(se) && length(se) > 0L) {
        if (is.null(names(se)) || any(!nzchar(names(se)))) {
          names(se) <- names(est)[seq_along(se)]
        }
        hit <- intersect(names(est), names(se))
        se_vec[hit] <- as.numeric(se[hit])
      }
      if (printSE) df$SE <- as.numeric(se_vec)
      if (print95CI) {
        df$`95% CI` <- ifelse(
          is.na(se_vec),
          NA_character_,
          paste0(
            "(", round(df$Estimate - stats::qnorm(0.975) * se_vec, 4),
            ", ", round(df$Estimate + stats::qnorm(0.975) * se_vec, 4), ")"
          )
        )
      }
      if (true_value_verbose && !is.null(true_values)) {
        tv <- .true_values_for_stage(true_values, stage)
        tv_vec <- rep(NA_real_, length(est)); names(tv_vec) <- names(est)
        if (!is.null(tv) && length(tv) > 0L) {
          if (is.null(names(tv)) || any(!nzchar(names(tv)))) {
            names(tv) <- names(est)[seq_along(tv)]
          }
          hit <- intersect(names(est), names(tv))
          tv_vec[hit] <- as.numeric(tv[hit])
        }
        df$`True value` <- as.numeric(tv_vec)
        df$Bias <- df$Estimate - df$`True value`
        df$RMSE <- sqrt(df$Bias^2)
      }
      rows[[length(rows) + 1L]] <- df
    }
    tables[[stage]] <- .bind_comparison_rows(rows, printSE, print95CI,
                                             true_value_verbose && !is.null(true_values))
  }

  out <- list(
    fits = fits,
    tables = tables,
    errors = errors,
    methods = methods,
    call = match.call()
  )
  class(out) <- "CIMEHR_method_comparisons"
  out
}

#' @export
print.CIMEHR_method_comparisons <- function(x, digits = 4, ...) {
  for (stage in c("visiting", "observation", "outcome")) {
    cat("\n---------------------------------\n")
    cat(stage, "process\n")
    cat("---------------------------------\n")
    print(round_numeric_df(x$tables[[stage]], digits = digits), row.names = FALSE)
  }
  if (length(x$errors) > 0L) {
    cat("\nMethods with errors:\n")
    print(x$errors)
  }
  invisible(x)
}

round_numeric_df <- function(x, digits = 4) {
  for (j in seq_along(x)) {
    if (is.numeric(x[[j]])) x[[j]] <- round(x[[j]], digits)
  }
  x
}

.empty_comparison_row <- function(method) {
  data.frame(
    Method = method,
    `Parameter of interest` = NA_character_,
    Estimate = NA_real_,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

.bind_comparison_rows <- function(rows, printSE, print95CI, include_truth) {
  cols <- c("Method", "Parameter of interest", "Estimate")
  if (printSE) cols <- c(cols, "SE")
  if (print95CI) cols <- c(cols, "95% CI")
  if (include_truth) cols <- c(cols, "True value", "Bias", "RMSE")

  rows <- lapply(rows, function(x) {
    missing <- setdiff(cols, names(x))
    for (nm in missing) x[[nm]] <- if (nm %in% c("Method", "Parameter of interest", "95% CI")) NA_character_ else NA_real_
    x[, cols, drop = FALSE]
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

.resolve_comparison_method <- function(method) {
  key <- tolower(gsub("[-. ]", "_", method))
  map <- list(
    cimehr = list(display = "CIMEHR", fun = "CIMEHR"),
    cimehr_timevarying_integral = list(display = "CIMEHR_timevarying_integral", fun = "CIMEHR_timevarying_integral"),
    cimehr_timevarying_ghq = list(display = "CIMEHR_timevarying_integral", fun = "CIMEHR_timevarying_integral"),
    cimehr_timevarying_ou = list(display = "CIMEHR_timevarying_ou", fun = "CIMEHR_timevarying_ou"),
    ehrjoint = list(display = "EHRJoint", fun = "EHRJoint"),
    inverse_intensity_rate_ratio_balancing = list(display = "Inverse_intensity_rate_ratio_balancing", fun = "Inverse_intensity_rate_ratio_balancing"),
    inverse_intensity_rate_ratio_weighting = list(display = "Inverse_intensity_rate_ratio_weighting", fun = "Inverse_intensity_rate_ratio_weighting"),
    joint_modeling_visiting_and_longitudinal_liang = list(display = "Joint_modeling_visiting_and_longitudinal_Liang", fun = "Joint_modeling_visiting_and_longitudinal_Liang"),
    linear_increment = list(display = "Linear_increment", fun = "Linear_increment_IP"),
    linear_increment_ip = list(display = "Linear_increment", fun = "Linear_increment_IP"),
    linear_mixed_model = list(display = "Linear_mixed_model", fun = "Linear_mixed_model"),
    multiple_imputation = list(display = "Multiple_imputation", fun = "Multiple_imputation_IP"),
    multiple_imputation_ip = list(display = "Multiple_imputation", fun = "Multiple_imputation_IP"),
    pairwise_likelihood = list(display = "Pairwise_likelihood", fun = "Pairwise_likelihood"),
    summary_stat = list(display = "Summary_stat", fun = "Summary_stat")
  )
  map[[key]]
}

.method_specific_user_args <- function(display_name, fun_name, method_args) {
  candidates <- unique(c(display_name, fun_name))
  for (nm in candidates) {
    if (!is.null(method_args[[nm]])) return(method_args[[nm]])
  }
  list()
}

.build_method_args <- function(display_name,
                               fun_name,
                               data,
                               extra_args,
                               method_args,
                               default_formula,
                               default_random_formula) {
  user_args <- .method_specific_user_args(display_name, fun_name, method_args)
  args_m <- list(data = data)
  if (length(extra_args) > 0L) {
    for (nm in names(extra_args)) args_m[[nm]] <- extra_args[[nm]]
  }
  if (length(user_args) > 0L) {
    for (nm in names(user_args)) args_m[[nm]] <- user_args[[nm]]
  }

  # Supply safe defaults for methods that otherwise require formulas. User
  # choices in method_args always override these defaults.
  if (identical(fun_name, "Linear_mixed_model")) {
    if (is.null(args_m$fixed)) args_m$fixed <- default_formula
    if (is.null(args_m$random)) args_m$random <- default_random_formula
    if (is.null(args_m$obs_indicator)) args_m$obs_indicator <- "R"
  }
  if (identical(fun_name, "Pairwise_likelihood")) {
    if (is.null(args_m$formula)) args_m$formula <- default_formula
    if (is.null(args_m$pair_sample)) args_m$pair_sample <- 5000L
  }
  if (identical(fun_name, "Summary_stat")) {
    # Summary_stat has usable defaults; no formula is needed.
  }
  if (identical(fun_name, "Linear_increment_IP")) {
    if (is.null(args_m$outcome_formula)) args_m$outcome_formula <- default_formula
    if (is.null(args_m$random_formula)) args_m$random_formula <- default_random_formula
  }
  if (identical(fun_name, "Multiple_imputation_IP")) {
    if (is.null(args_m$outcome_formula)) args_m$outcome_formula <- default_formula
    if (is.null(args_m$random_formula)) args_m$random_formula <- default_random_formula
  }

  args_m
}

# User-facing aliases requested for the comparison workflow.
#' @rdname method_comparisons
#' @export
Linear_increment <- function(...) Linear_increment_IP(...)

#' @rdname method_comparisons
#' @export
Multiple_imputation <- function(...) Multiple_imputation_IP(...)


#' Fit and Compare Multiple Methods
#'
#' Fits any combination of implemented CIMEHR package methods and returns
#' process-specific comparison tables. Each table contains one row per method
#' and parameter. For methods that do not estimate a requested process, the
#' table contains an \code{NA} row for that process.
#'
#' @param data Long-format analysis data. This can also be a \code{sim_CIMEHR}
#'   object returned by \code{sim_data_gen()}, in which case \code{data$long_data}
#'   is used and true simulation parameters are detected automatically.
#' @param methods Character vector of method function names to fit.
#' @param method_args Optional named list of argument lists. Names should match
#'   entries in \code{methods}. Each sub-list is passed to the corresponding
#'   fitting function.
#' @param printSE Logical; include standard error column where available.
#' @param print95CI Logical; include Wald 95 percent confidence interval columns
#'   where standard errors are available.
#' @param true_values Optional named list of true parameters. If \code{NULL},
#'   the function tries to use \code{attr(data, "true_params")} or
#'   \code{data$params} for simulated data.
#' @param true_value_verbose Logical; include true value and bias columns when
#'   true values are available.
#' @param which_balancing Character; for
#'   \code{Inverse_intensity_rate_ratio_balancing}, choose outcome estimates
#'   from \code{"balancing"} or \code{"mle"} weights.
#' @param stop_on_error Logical; if \code{TRUE}, stop when a method fails. If
#'   \code{FALSE}, store the error and continue fitting the remaining methods.
#' @param ... Additional arguments passed to all methods unless overridden by
#'   method-specific entries in \code{method_args}.
#'
#' @return An object of class \code{"CIMEHR_method_comparisons"} with fitted
#'   model objects and process-specific tables.
#' @export
method_comparisons <- function(data,
                               methods = c("CIMEHR", "CIMEHR_timevarying_integral", "CIMEHR_timevarying_ou"),
                               method_args = list(),
                               printSE = TRUE,
                               print95CI = TRUE,
                               true_values = NULL,
                               true_value_verbose = TRUE,
                               which_balancing = c("balancing", "mle"),
                               stop_on_error = FALSE,
                               ...) {
  which_balancing <- match.arg(which_balancing)
  extra_args <- list(...)

  # Allow users to pass the full object returned by sim_data_gen().
  if (inherits(data, "sim_CIMEHR") && !is.null(data$long_data)) {
    if (is.null(true_values)) true_values <- data$params
    data <- data$long_data
  }
  if (is.null(true_values)) true_values <- attr(data, "true_params")

  fits <- list()
  errors <- list()

  for (m in methods) {
    if (!exists(m, mode = "function")) {
      msg <- paste0("Method function not found: ", m)
      if (stop_on_error) stop(msg, call. = FALSE)
      warning(msg, call. = FALSE)
      fits[[m]] <- NULL
      errors[[m]] <- msg
      next
    }
    f <- get(m, mode = "function")
    args_m <- c(list(data = data), extra_args, method_args[[m]])
    fit_m <- tryCatch(do.call(f, args_m), error = function(e) e)
    if (inherits(fit_m, "error")) {
      if (stop_on_error) stop(fit_m)
      warning("Method ", m, " failed: ", conditionMessage(fit_m), call. = FALSE)
      fits[[m]] <- NULL
      errors[[m]] <- conditionMessage(fit_m)
    } else {
      fits[[m]] <- fit_m
    }
  }

  stages <- c("visiting", "observation", "outcome")
  tables <- setNames(vector("list", length(stages)), stages)
  for (stage in stages) {
    rows <- list()
    for (m in methods) {
      fit <- fits[[m]]
      if (is.null(fit)) {
        rows[[length(rows) + 1L]] <- data.frame(Method = m, Parameter = NA_character_, Estimate = NA_real_, stringsAsFactors = FALSE)
        next
      }
      if (!.method_has_stage(fit, stage)) {
        rows[[length(rows) + 1L]] <- data.frame(Method = m, Parameter = NA_character_, Estimate = NA_real_, stringsAsFactors = FALSE)
        next
      }
      ex <- .extract_stage_estimates(fit, stage = stage, which_estimate = which_balancing)
      est <- ex$estimate
      se <- ex$se
      if (length(est) == 0L) {
        rows[[length(rows) + 1L]] <- data.frame(Method = m, Parameter = NA_character_, Estimate = NA_real_, stringsAsFactors = FALSE)
        next
      }
      if (is.null(names(est)) || any(names(est) == "")) names(est) <- paste0("param", seq_along(est))
      df <- data.frame(Method = m, Parameter = names(est), Estimate = as.numeric(est), stringsAsFactors = FALSE)
      if (printSE) {
        se_vec <- rep(NA_real_, length(est)); names(se_vec) <- names(est)
        if (!is.null(se)) {
          if (is.null(names(se))) names(se) <- names(est)[seq_along(se)]
          hit <- intersect(names(est), names(se))
          se_vec[hit] <- as.numeric(se[hit])
        }
        df$SE <- as.numeric(se_vec)
      }
      if (print95CI) {
        if (!"SE" %in% names(df)) {
          se_vec <- rep(NA_real_, length(est)); names(se_vec) <- names(est)
          if (!is.null(se)) {
            if (is.null(names(se))) names(se) <- names(est)[seq_along(se)]
            hit <- intersect(names(est), names(se))
            se_vec[hit] <- as.numeric(se[hit])
          }
        } else {
          se_vec <- df$SE
        }
        df$CI_low <- df$Estimate - stats::qnorm(0.975) * se_vec
        df$CI_high <- df$Estimate + stats::qnorm(0.975) * se_vec
      }
      if (true_value_verbose && !is.null(true_values)) {
        tv <- .true_values_for_stage(true_values, stage)
        tv_vec <- rep(NA_real_, length(est)); names(tv_vec) <- names(est)
        if (!is.null(tv)) {
          if (is.null(names(tv))) names(tv) <- names(est)[seq_along(tv)]
          hit <- intersect(names(est), names(tv))
          tv_vec[hit] <- as.numeric(tv[hit])
        }
        df$True_value <- as.numeric(tv_vec)
        df$Bias <- df$Estimate - df$True_value
      }
      rows[[length(rows) + 1L]] <- df
    }
    tables[[stage]] <- do.call(rbind, rows)
    rownames(tables[[stage]]) <- NULL
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

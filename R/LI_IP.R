#' Linear Increment Imputation with Informative Presence Adjustment
#'
#' Combines linear increment (LI) imputation for missing outcomes with IP-adjusted
#' analysis methods. First imputes missing Y values using a linear increment model
#' (Aalen & Gunnes, 2010), then applies an IP-adjustment method to the completed dataset.
#'
#' @param data Data frame in long format with one row per visit.
#' @param id Character. Subject ID column name. Default is "id".
#' @param time Character. Time variable column name. Default is "time".
#' @param outcome Character. Outcome variable column name. Default is "Y".
#' @param obs_indicator Character. Observation indicator column name
#'        (1 = outcome observed, 0 = outcome missing at visit). Default is "R".
#'        Rows with NA in time are treated as unvisited time points.
#'
#' @section Linear Increment Model:
#' The linear increment model (Aalen & Gunnes, 2010) models the change in outcome
#' between consecutive observations rather than the outcome level itself. This approach
#' is particularly suited for irregular longitudinal data because:
#' \itemize{
#'   \item It naturally handles varying time gaps between observations
#'   \item It assumes the increment \eqn{\Delta Y = Y(t_j) - Y(t_{j-1})} follows a
#'         linear model with variance proportional to the time gap
#'   \item Missing values are imputed by predicting the expected increment
#' }
#'
#' The model can be written as:
#' \deqn{\Delta Y_{ij} = (t_j - t_{j-1}) \cdot (\beta_0 + \beta' X_{ij}) + \epsilon_{ij}}
#' where \eqn{Var(\epsilon_{ij}) \propto (t_j - t_{j-1})}.
#'
#' @param li_formula Formula. Model formula for the linear increment model.
#'        Default is Y ~ Z + X + time. The left-hand side should be the outcome variable.
#' @param li_covariance Character. Covariance structure for the linear increment model.
#'        One of:
#'   \itemize{
#'     \item "randomwalk" - Random walk covariance (default, recommended for irregular data)
#'     \item "independence" - Independent increments
#'     \item "exchangeable" - Exchangeable correlation
#'   }
#' @param li_limit Formula. Limit specification for the linear increment model.
#'        Default is ~ 1 (constant limit). See \code{\link[slim]{slim}} for details.
#'
#' @section IP Method Parameters:
#' @param ip_method Character. IP adjustment method. One of:
#'   \itemize{
#'     \item "LMM" - Standard linear mixed model (no IP adjustment, for comparison)
#'     \item "Liang" - Joint model of Liang et al. (2009)
#'     \item "IIRR_weighting" - Inverse intensity rate ratio weighting
#'     \item "IIRR_balancing" - IIRR with balancing weights (Yiu & Su, 2025)
#'     \item "PairLik" - Pairwise composite likelihood
#'   }
#' @param visit_covs Character vector. Covariates for visiting process model.
#'        Used by Liang, IIRR methods. Default is c("Z", "X").
#' @param outcome_formula Formula. Model formula for longitudinal outcome.
#'        Used by LMM, IIRR_balancing, PairLik. Default is Y ~ Z + X + time.
#' @param random_formula Formula. Random effects formula for LMM.
#'        Default is ~ 1 | id.
#' @param outcome_covs Character vector. Covariates for outcome model in Liang method.
#'        Default is NULL (uses visit_covs).
#' @param random_covs Character vector. Covariates with random effects in Liang method.
#'        Default is first element of visit_covs.
#' @param ip_args List. Additional arguments passed to the IP method. See
#'        \code{\link{MI_IP_fun}} for detailed options for each method. Common options:
#'   \describe{
#'     \item{IIRR methods:}{\code{phi} - sensitivity parameter (default 0)}
#'     \item{PairLik:}{\code{pair_sample} - number of pairs to sample}
#'     \item{LMM:}{\code{visit_adjust} - "OA" or "VA" for visit adjustment}
#'   }
#'
#' @return A list containing:
#'   \item{coefficients}{Named vector of coefficient estimates.}
#'   \item{ip_method}{The IP method used.}
#'   \item{li_model}{Summary of the linear increment model used for imputation.}
#'   \item{n_imputed}{Number of values imputed.}
#'   \item{imputed_data}{The completed dataset after LI imputation.}
#'   \item{call}{The function call.}
#'
#' @details
#' Unlike multiple imputation (MI), the linear increment approach produces a single
#' completed dataset. This has implications:
#' \itemize{
#'   \item No pooling of estimates is needed
#'   \item Standard errors from the IP method do not account for imputation uncertainty
#'   \item Computationally faster than MI
#' }
#'
#' The LI method is appropriate when:
#' \itemize
#'   \item The outcome evolves smoothly over time
#'   \item Missing data is MAR (missing at random) given observed history
#'   \item Time gaps between observations vary across subjects
#' }
#'
#' @references
#' Aalen, O. O., & Gunnes, N. (
#' 2010). A dynamic approach for reconstructing missing
#' longitudinal data using the linear increments model. Biostatistics, 11(3), 453-472.
#'
#' Liang, Y., Lu, W., & Ying, Z. (2009). Joint modeling and analysis of
#' longitudinal data with informative observation times. Biometrics, 65, 377-384.
#'
#' @examples
#' \dontrun{
#' # Example 1: Liang's method with default LI imputation
#' result <- LI_IP(
#'   data = mydata,
#'   ip_method = "Liang",
#'   visit_covs = c("Z", "X")
#' )
#'
#' # Example 2: Custom LI model formula
#' result <- LI_IP(
#'   data = mydata,
#'   ip_method = "IIRR_weighting",
#'   li_formula = Y ~ Z + X + time + Z:time,
#'   li_covariance = "randomwalk",
#'   visit_covs = c("Z", "X")
#' )
#'
#' # Example 3: IIRR with sensitivity analysis (phi != 0)
#' result <- LI_IP(
#'   data = mydata,
#'   ip_method = "IIRR_balancing",
#'   visit_covs = c("Z", "X"),
#'   ip_args = list(phi = 0.3)  # Sensitivity parameter
#' )
#'
#' # Example 4: Compare LI vs MI approaches
#' result_li <- LI_IP(
#'   data = mydata,
#'   ip_method = "Liang",
#'   visit_covs = c("Z", "X")
#' )
#' result_mi <- MI_IP_fun(
#'   data = mydata,
#'   ip_method = "Liang",
#'   visit_covs = c("Z", "X"),
#'   m = 5
#' )
#'
#' # Compare coefficients
#' rbind(LI = result_li$coefficients, MI = result_mi$coefficients)
#' }
#'
#' @seealso \code{\link{MI_IP_fun}} for multiple imputation approach
#'
#' @export
LI_IP <- function(
    data,
    # Column names
    id = "id",
    time = "time",
    outcome = "Y",
    obs_indicator = "R",

    # Linear Increment model settings
    li_formula = NULL,
    li_covariance = c("randomwalk", "independence", "exchangeable"),
    li_limit = ~ 1,

    # IP method settings
    ip_method = c("LMM", "Liang", "IIRR_weighting", "IIRR_balancing", "PairLik"),

    # Common covariate specifications
    visit_covs = c("Z", "X"),
    outcome_covs = NULL,
    random_covs = NULL,

    # Formula specifications (for IP methods that use them)
    outcome_formula = NULL,
    random_formula = NULL,

    # Method-specific arguments
    ip_args = list()
) {

  # ============================================
  # Input validation and setup
  # ============================================
  call <- match.call()
  li_covariance <- match.arg(li_covariance)
  ip_method <- match.arg(ip_method)

  # Check required packages
  if (!requireNamespace("slim", quietly = TRUE)) {
    stop("Package 'slim' is required. Install with: install.packages('slim')")
  }
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Package 'data.table' is required. Install with: install.packages('data.table')")
  }

  # Set defaults
  if (is.null(outcome_covs)) outcome_covs <- visit_covs
  if (is.null(random_covs)) random_covs <- visit_covs[1]

  # Default formulas
  if (is.null(li_formula)) {
    li_formula <- as.formula(
      paste(outcome, "~", paste(c(visit_covs, time), collapse = " + "))
    )
  }
  if (is.null(outcome_formula)) {
    outcome_formula <- as.formula(
      paste(outcome, "~", paste(c(outcome_covs, time), collapse = " + "))
    )
  }
  if (is.null(random_formula)) {
    random_formula <- as.formula(paste("~ 1 |", id))
  }

  # ============================================
  # Stage 1: Linear Increment Imputation
  # ============================================
  message("Stage 1: Linear Increment Imputation...")
  message("  LI model: ", deparse(li_formula))
  message("  Covariance structure: ", li_covariance)

  li_result <- .li_impute_outcome(
    data = data,
    id = id,
    time = time,
    outcome = outcome,
    obs_indicator = obs_indicator,
    formula = li_formula,
    covariance = li_covariance,
    limit = li_limit
  )

  imputed_data <- li_result$data
  n_imputed <- li_result$n_imputed

  message("  Imputed ", n_imputed, " missing values")

  # ============================================
  # Stage 2: Apply IP method
  # ============================================
  message("\nStage 2: Applying ", ip_method, "...")

  fit <- tryCatch({
    .apply_ip_method(
      data = imputed_data,
      ip_method = ip_method,
      id = id,
      time = time,
      outcome = outcome,
      obs_indicator = obs_indicator,
      visit_covs = visit_covs,
      outcome_covs = outcome_covs,
      random_covs = random_covs,
      outcome_formula = outcome_formula,
      random_formula = random_formula,
      ip_args = ip_args
    )
  }, error = function(e) {
    stop("IP method failed: ", e$message)
  })

  # ============================================
  # Return results
  # ============================================
  message("Done!")

  structure(
    list(
      coefficients = fit$coef,
      se = fit$se,
      ip_method = ip_method,
      ip_args = ip_args,
      li_model = list(
        formula = li_formula,
        covariance = li_covariance,
        limit = li_limit
      ),
      n_imputed = n_imputed,
      imputed_data = imputed_data,
      call = call
    ),
    class = "LI_IP_result"
  )
}


# ============================================
# Internal: Linear Increment Imputation
# ============================================
#' @description
#' Linear Increment Model for imputation:
#' The model assumes that increments Delta_Y = Y(t_j) - Y(t_{j-1}) follow:
#'   Delta_Y_ij = (t_j - t_{j-1}) * (beta_0 + beta' X_ij) + epsilon_ij
#' where Var(epsilon_ij) is proportional to (t_j - t_{j-1}).
#'
#' This "random walk" structure means larger time gaps have more uncertainty,
#' which is natural for irregular longitudinal data.
#' @keywords internal
.li_impute_outcome <- function(data, id, time, outcome, obs_indicator,
                               formula, covariance, limit) {

  # slim::slim requires a keyed data.table
  # Convert to data.table and set key
  dt <- data.table::as.data.table(data)

  # Set key using the id and time columns (must use setkeyv for character column names)
  data.table::setkeyv(dt, c(id, time))

  # Identify rows for fitting: need time and all covariates and observed Y
  formula_vars <- all.vars(formula)

  fit_idx <- !is.na(dt[[time]]) & !is.na(dt[[outcome]])
  for (v in formula_vars) {
    if (v != outcome && v %in% colnames(dt)) {
      fit_idx <- fit_idx & !is.na(dt[[v]])
    }
  }

  # Create keyed data.table for fitting
  dt_fit <- dt[fit_idx]
  data.table::setkeyv(dt_fit, c(id, time))

  if (nrow(dt_fit) < 10) {
    stop("Too few complete observations to fit LI model (n = ", nrow(dt_fit), ")")
  }

  # Fit SLIM (linear increment model) - requires keyed data.table
  fit_slim <- slim::slim(
    formula = formula,
    data = dt_fit,
    covariance = covariance,
    limit = limit
  )

  # Identify rows to predict: time is known but Y is missing
  pred_idx <- !is.na(dt[[time]]) & is.na(dt[[outcome]])
  for (v in formula_vars) {
    if (v != outcome && v %in% colnames(dt)) {
      pred_idx <- pred_idx & !is.na(dt[[v]])
    }
  }

  n_imputed <- sum(pred_idx)

  # Predict missing Y values and fill in
  if (n_imputed > 0) {
    # Create keyed data.table for prediction
    dt_pred <- dt[pred_idx]
    data.table::setkeyv(dt_pred, c(id, time))

    # Get predictions for missing values
    pred_values <- predict(fit_slim, newdata = dt_pred)

    # Fill in predictions where outcome is missing
    dt[[outcome]][pred_idx] <- pred_values
  }

  list(
    data = as.data.frame(dt),
    n_imputed = n_imputed,
    slim_fit = fit_slim
  )
}


# ============================================
# Internal: Apply IP Method (shared with MI_IP_fun)
# ============================================
# Note: This function is the same as in MI_IP_fun.R
# In a package, you would define it once in a utils.R file
.apply_ip_method <- function(data, ip_method, id, time, outcome, obs_indicator,
                             visit_covs, outcome_covs, random_covs,
                             outcome_formula, random_formula,
                             ip_args = list()) {

  result <- switch(ip_method,
                   "LMM" = {
                     # Extract LMM-specific arguments with defaults
                     lmm_method <- ip_args$method %||% "REML"
                     lmm_optimizer <- ip_args$optimizer %||% "nlminb"
                     lmm_visit_adjust <- ip_args$visit_adjust  # NULL is valid default

                     coefs <- LMM_fun(
                       data = data,
                       id = id,
                       fixed = outcome_formula,
                       random = random_formula,
                       obs_indicator = obs_indicator,
                       visit_adjust = lmm_visit_adjust,
                       method = lmm_method,
                       optimizer = lmm_optimizer
                     )
                     list(
                       coef = coefs[, "Value"],
                       se = coefs[, "Std.Error"]
                     )
                   },

                   "Liang" = {
                     fit <- JMVL_Liang_fun(
                       data = data,
                       id = id,
                       time = time,
                       outcome = outcome,
                       obs_indicator = obs_indicator,
                       visit_covs = visit_covs,
                       long_covs = outcome_covs,
                       random_covs = random_covs
                     )
                     list(
                       coef = fit$beta_hat,
                       se = NULL,  # Liang doesn't return SE directly
                       gamma_hat = fit$gamma_hat
                     )
                   },

                   "IIRR_weighting" = {
                     # Extract IIRR-specific arguments with defaults
                     phi <- ip_args$phi %||% 0
                     balance_covs <- ip_args$balance_covs  # NULL uses visit_covs in function

                     fit <- IIRR_balancing_fun(
                       data = data,
                       id = id,
                       time = time,
                       outcome = outcome,
                       obs_indicator = obs_indicator,
                       visit_covs = visit_covs,
                       balance_covs = balance_covs,
                       formula = outcome_formula,
                       phi = phi
                     )
                     # Use MLE weights (standard IIRR)
                     list(
                       coef = fit$beta_mle,
                       se = NULL,
                       gamma_hat = fit$gamma_hat,
                       phi = phi
                     )
                   },

                   "IIRR_balancing" = {
                     # Extract IIRR-specific arguments with defaults
                     phi <- ip_args$phi %||% 0
                     balance_covs <- ip_args$balance_covs  # NULL uses visit_covs in function

                     fit <- IIRR_balancing_fun(
                       data = data,
                       id = id,
                       time = time,
                       outcome = outcome,
                       obs_indicator = obs_indicator,
                       visit_covs = visit_covs,
                       balance_covs = balance_covs,
                       formula = outcome_formula,
                       phi = phi
                     )
                     # Use balancing weights
                     list(
                       coef = fit$beta_bal,
                       se = NULL,
                       gamma_hat = fit$gamma_hat,
                       phi = phi
                     )
                   },

                   "PairLik" = {
                     # Extract PairLik-specific arguments with defaults
                     pair_sample <- ip_args$pair_sample  # NULL uses all pairs
                     pair_family <- ip_args$family %||% "gaussian"

                     fit <- pairwise_glm_fit_rcpp(
                       data = data[data[[obs_indicator]] == 1, ],
                       formula = outcome_formula,
                       id_col = id,
                       time_col = time,
                       y_col = outcome,
                       family = pair_family,
                       pair_sample = pair_sample
                     )
                     list(
                       coef = fit$coef,
                       se = if (!is.null(fit$vcov_naive)) sqrt(diag(fit$vcov_naive)) else NULL
                     )
                   },

                   stop("Unknown IP method: ", ip_method)
  )

  result
}


# ============================================
# Print method
# ============================================
#' @export
print.LI_IP_result <- function(x, ...) {
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


# ============================================
# Summary method
# ============================================
#' @export
summary.LI_IP_result <- function(object, ...) {
  cat("Linear Increment Imputation with IP Adjustment\n")
  cat("===============================================\n\n")

  cat("Call:\n")
  print(object$call)
  cat("\n")

  cat("Linear Increment Model:\n")
  cat("  Formula:", deparse(object$li_model$formula), "\n")
  cat("  Covariance structure:", object$li_model$covariance, "\n")
  cat("  Number of values imputed:", object$n_imputed, "\n\n")

  cat("IP Method:", object$ip_method, "\n\n")

  cat("Coefficients:\n")
  if (!is.null(object$se)) {
    z_val <- object$coefficients / object$se
    p_val <- 2 * pnorm(-abs(z_val))
    print_df <- data.frame(
      Estimate = object$coefficients,
      `Std.Error` = object$se,
      `z value` = z_val,
      `Pr(>|z|)` = p_val,
      row.names = names(object$coefficients),
      check.names = FALSE
    )
  } else {
    print_df <- data.frame(
      Estimate = object$coefficients,
      row.names = names(object$coefficients)
    )
    cat("(Standard errors not available for this IP method)\n")
  }
  print(round(print_df, 4))

  invisible(object)
}

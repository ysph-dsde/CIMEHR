#' Multiple Imputation with Informative Presence Adjustment
#'
#' Combines multiple imputation for missing outcomes with IP-adjusted analysis
#' methods. First imputes missing Y values at observed visits, then applies
#' an IP-adjustment method to each imputed dataset and pools results using
#' Rubin's rules.
#'
#' @param data Data frame in long format with one row per visit.
#' @param id Character. Subject ID column name. Default is "id".
#' @param time Character. Time variable column name. Default is "time".
#' @param outcome Character. Outcome variable column name. Default is "Y".
#' @param obs_indicator Character. Observation indicator column name
#'        (1 = outcome observed, 0 = outcome missing at visit). Default is "R".
#'        Rows with NA in time are treated as unvisited time points.
#'
#' @section Imputation Parameters:
#' @param m Integer. Number of imputations. Default is 5.
#' @param maxit Integer. Maximum iterations for mice. Default is 50.
#' @param impute_method Character. Method for imputing Y. One of:
#'   \itemize{
#'     \item "2l.pan" - Two-level model with heterogeneous within-group variance (default)
#'     \item "2l.lmer" - Two-level model using lmer
#'     \item "2l.norm" - Two-level normal model
#'   }
#' @param impute_fixed Character vector. Fixed effect covariates for imputation model.
#'        Default is NULL, which uses c(visit_covs, time).
#'        The imputation model for Y is:
#'        \deqn{Y_{ij} = \beta_0 + \sum_k \beta_k X_{kij} + b_{0i} + \sum_l b_{li} Z_{lij} + \epsilon_{ij}}
#'        where fixed effects include an intercept plus all variables in impute_fixed.
#' @param impute_random Character vector. Random slope variables for imputation model.
#'        Default is NULL (random intercept only). For example, c("Z") adds a random
#'        slope for Z: \eqn{b_{0i} + b_{1i} Z_{ij}}.
#' @param impute_include_R Logical. Whether to include observation indicator R as a
#'        predictor in the imputation model. Default is TRUE. This is recommended
#'        because R may be informative about missing Y values (MAR assumption).
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
#' @param ip_args List. Additional arguments passed to the IP method. Method-specific
#'        options include:
#'   \describe{
#'     \item{IIRR_weighting, IIRR_balancing:}{
#'       \itemize{
#'         \item \code{phi}: Sensitivity parameter (default 0). phi = 0 assumes
#'               Visiting at Random (VAR); phi > 0 means higher outcome increases
#'               visit probability; phi < 0 means higher outcome decreases visit probability.
#'         \item \code{balance_covs}: Character vector of covariates to balance
#'               (default uses visit_covs).
#'       }
#'     }
#'     \item{PairLik:}{
#'       \itemize{
#'         \item \code{pair_sample}: Integer. Number of subject pairs to sample
#'               for large datasets (default NULL uses all pairs).
#'         \item \code{family}: "gaussian" (default) or "binomial".
#'       }
#'     }
#'     \item{LMM:}{
#'       \itemize{
#'         \item \code{method}: "REML" (default) or "ML".
#'         \item \code{optimizer}: "nlminb" (default) or "optim".
#'         \item \code{visit_adjust}: NULL, "OA", or "VA" for visit count adjustment.
#'       }
#'     }
#'   }
#'
#' @section Pooling Parameters:
#' @param pool_method Character. Method for pooling estimates. Currently only
#'        "mean" is implemented (simple averaging of point estimates).
#'        Future versions will support "rubin" for proper Rubin's rules with SE.
#'
#' @return A list containing:
#'   \item{coefficients}{Named vector of pooled coefficient estimates.}
#'   \item{by_imputation}{List of results from each imputed dataset.}
#'   \item{ip_method}{The IP method used.}
#'   \item{m}{Number of imputations.}
#'   \item{call}{The function call.}
#'
#' @details
#' The function implements a two-stage procedure:
#'
#' \strong{Stage 1: Multiple Imputation}
#'
#' Missing outcome values (where R = 0 but time is not NA) are imputed using
#' a multilevel model that accounts for the clustered structure of longitudinal
#' data. The imputation model is:
#' \deqn{Y_{ij} = X_{ij}'\beta + Z_{ij}'b_i + \epsilon_{ij}}
#' where \eqn{b_i \sim N(0, \Sigma_b)} and \eqn{\epsilon_{ij} \sim N(0, \sigma^2)}.
#'
#' \strong{Stage 2: IP-Adjusted Analysis}
#'
#' Each imputed dataset is analyzed using the specified IP method:
#'
#' \emph{Liang et al. (2009):} Joint models the visiting process and outcome:
#' \deqn{\lambda_i(t) = \eta_i \exp(\gamma' X_i^V) \lambda_0(t)}
#' \deqn{Y_i(t) = \beta' X_i^Y(t) + b_i' Z_i^Y(t) + \epsilon_i(t)}
#' with \eqn{E[b_i | \eta_i] = \theta(\eta_i - 1)} linking the processes.
#'
#' \emph{IIRR Weighting:} Weights observations by inverse visit intensity:
#' \deqn{w_i(t) = 1/\hat{\lambda}_i(t)}
#'
#' \emph{IIRR Balancing:} Solves for weights satisfying balance constraints:
#' \deqn{\sum_i w_i h_i(t) = \sum_i h_i(t) d\hat{\Lambda}_0(t)}
#'
#' @references
#' Liang, Y., Lu, W., & Ying, Z. (2009). Joint modeling and analysis of
#' longitudinal data with informative observation times. Biometrics, 65, 377-384.
#'
#' Yiu, S. & Su, L. (2025). Accommodating informative visit times: a sensitivity
#' analysis approach with balancing weights estimators. JRSS-C, 74, 824-843.
#'
#' Rubin, D. B. (1987). Multiple Imputation for Nonresponse in Surveys. Wiley.
#'
#' @examples
#' \dontrun{
#' # Example 1: Liang's method with default imputation
#' # Imputation model: Y ~ Z + X + time + (1 | id)
#' result <- MI_IP(
#'   data = mydata,
#'   ip_method = "Liang",
#'   visit_covs = c("Z", "X"),
#'   m = 5
#' )
#'
#' # Example 2: Specify imputation model explicitly
#' # Imputation model: Y ~ Z + X + time + R + (1 + Z | id)
#' result <- MI_IP(
#'   data = mydata,
#'   ip_method = "Liang",
#'   visit_covs = c("Z", "X"),
#'   impute_fixed = c("Z", "X", "time"),
#'   impute_random = "Z",
#'   impute_include_R = TRUE,
#'   m = 5
#' )
#'
#' # Example 3: IIRR balancing with sensitivity analysis (phi != 0)
#' result <- MI_IP(
#'   data = mydata,
#'   ip_method = "IIRR_balancing",
#'   visit_covs = c("Z", "X"),
#'   outcome_formula = Y ~ Z + X + time,
#'   ip_args = list(phi = 0.3),  # Sensitivity parameter
#'   m = 10
#' )
#'
#' # Example 4: Standard LMM with visit adjustment
#' result <- MI_IP(
#'   data = mydata,
#'   ip_method = "LMM",
#'   outcome_formula = Y ~ Z + X + time,
#'   random_formula = ~ 1 + Z | id,
#'   ip_args = list(
#'     method = "REML",
#'     visit_adjust = "VA"  # Visit-aligned adjustment
#'   )
#' )
#'
#' # Example 5: Pairwise likelihood with subsampling for large data
#' result <- MI_IP(
#'   data = large_data,
#'   ip_method = "PairLik",
#'   outcome_formula = Y ~ Z + X + time,
#'   ip_args = list(pair_sample = 5000),  # Sample 5000 pairs
#'   m = 5
#' )
#' }
#'
#' @export
MI_IP <- function(
    data,
    # Column names
    id = "id",
    time = "time",
    outcome = "Y",
    obs_indicator = "R",

    # Imputation settings
    m = 5,
    maxit = 50,
    impute_method = c("2l.pan", "2l.lmer", "2l.norm"),
    impute_fixed = NULL,
    impute_random = NULL,
    impute_include_R = TRUE,


    # IP method settings
    ip_method = c("LMM", "Liang", "IIRR_weighting", "IIRR_balancing", "PairLik"),

    # Common covariate specifications
    visit_covs = c("Z", "X"),
    outcome_covs = NULL,
    random_covs = NULL,

    # Formula specifications (for methods that use them)
    outcome_formula = NULL,
    random_formula = NULL,

    # Method-specific arguments
    ip_args = list(),

    # Pooling
    pool_method = c("mean", "rubin")
) {

  # ============================================
  # Input validation and setup
  # ============================================
  call <- match.call()
  impute_method <- match.arg(impute_method)
  ip_method <- match.arg(ip_method)
  pool_method <- match.arg(pool_method)

  # Check required packages
  if (!requireNamespace("mice", quietly = TRUE)) {
    stop("Package 'mice' is required. Install with: install.packages('mice')")
  }
  if (!requireNamespace("nlme", quietly = TRUE)) {
    stop("Package 'nlme' is required. Install with: install.packages('nlme')")
  }

  # Set defaults
  if (is.null(outcome_covs)) outcome_covs <- visit_covs
  if (is.null(random_covs)) random_covs <- visit_covs[1]
  if (is.null(impute_fixed)) {
    impute_fixed <- unique(c(visit_covs, time))
  }

  # Default formulas
  if (is.null(outcome_formula)) {
    outcome_formula <- as.formula(
      paste(outcome, "~", paste(c(outcome_covs, time), collapse = " + "))
    )
  }
  if (is.null(random_formula)) {
    random_formula <- as.formula(paste("~ 1 |", id))
  }

  # ============================================
  # Stage 1: Multiple Imputation
  # ============================================
  message("Stage 1: Imputing missing outcomes...")
  message("  Imputation model: ", outcome, " ~ ",
          paste(impute_fixed, collapse = " + "),
          " + (1",
          if (!is.null(impute_random)) paste0(" + ", paste(impute_random, collapse = " + ")),
          " | ", id, ")")

  imputed_data <- .mi_impute_outcome(
    data = data,
    id = id,
    time = time,
    outcome = outcome,
    obs_indicator = obs_indicator,
    m = m,
    maxit = maxit,
    method = impute_method,
    fixed_vars = impute_fixed,
    random_vars = impute_random,
    include_R = impute_include_R
  )

  # ============================================
  # Stage 2: Apply IP method to each imputed dataset
  # ============================================
  message("\nStage 2: Applying ", ip_method, " to each imputed dataset...")

  fits <- lapply(seq_along(imputed_data), function(i) {
    message("  Processing imputation ", i, "/", length(imputed_data))
    tryCatch({
      .apply_ip_method(
        data = imputed_data[[i]],
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
      warning("IP method failed for imputation ", i, ": ", e$message)
      NULL
    })
  })

  # Remove failed fits
  fits <- Filter(Negate(is.null), fits)
  if (length(fits) == 0) {
    stop("All imputations failed. Check data and model specifications.")
  }

  # ============================================
  # Stage 3: Pool results
  # ============================================
  message("\nStage 3: Pooling results using ", pool_method, " method...")
  pooled <- .pool_estimates(fits, method = pool_method)

  # ============================================
  # Return results
  # ============================================
  message("Done!")

  structure(
    list(
      coefficients = pooled$coef,
      se = pooled$se,
      by_imputation = fits,
      ip_method = ip_method,
      ip_args = ip_args,
      impute_method = impute_method,
      impute_model = list(
        fixed = impute_fixed,
        random = impute_random,
        include_R = impute_include_R
      ),
      m = length(fits),
      m_requested = m,
      call = call
    ),
    class = "MI_IP_result"
  )
}


# ============================================
# Internal: Multiple Imputation
# ============================================
#' @description
#' Imputation model is a two-level mixed model:
#' Y_ij = beta_0 + sum(beta_k * X_kij) + b_0i + sum(b_li * Z_lij) + epsilon_ij
#' where:
#'   - X_kij are the fixed effect covariates (fixed_vars)
#'   - Z_lij are the random slope covariates (random_vars)
#'   - b_0i is the random intercept for subject i
#'   - b_li are the random slopes for subject i
#' @keywords internal
.mi_impute_outcome <- function(data, id, time, outcome, obs_indicator,
                               m, maxit, method, fixed_vars, random_vars,
                               include_R = TRUE) {

  data <- as.data.frame(data)

  # Separate unvisited rows (time is NA)
  unvisited <- data[is.na(data[[time]]), ]
  visited <- data[!is.na(data[[time]]), ]

  if (nrow(visited) == 0) {
    stop("No visited time points found (all time values are NA)")
  }

  # Build predictor list
  predictors <- fixed_vars
  if (include_R && obs_indicator %in% colnames(visited)) {
    predictors <- unique(c(predictors, obs_indicator))
  }

  # Setup mice
  imp0 <- mice::mice(visited, maxit = 0, print = FALSE)
  pred_matrix <- imp0$predictorMatrix
  meth <- imp0$method

  # Only impute outcome
  meth[] <- ""
  meth[outcome] <- method

  # Setup predictor matrix for outcome imputation
  # -2 indicates the cluster/grouping variable
  # 1 indicates a fixed effect predictor
  pred_matrix[, ] <- 0
  pred_matrix[outcome, id] <- -2  # Cluster variable

  for (v in predictors) {
    if (v %in% colnames(pred_matrix) && v != outcome && v != id) {
      pred_matrix[outcome, v] <- 1
    }
  }

  # Where to impute: Y is NA and time is not NA (i.e., visited but not observed)
  where_mat <- matrix(FALSE, nrow = nrow(visited), ncol = ncol(visited))
  colnames(where_mat) <- colnames(visited)
  where_mat[is.na(visited[[outcome]]) & !is.na(visited[[time]]), outcome] <- TRUE

  # Random slope specification via blots
  # This tells mice which variables should have random slopes in the imputation model
  blots <- list()
  if (!is.null(random_vars) && length(random_vars) > 0) {
    blots[[outcome]] <- list(random.slope = random_vars)
  }

  # Print model summary for user
  n_to_impute <- sum(where_mat[, outcome])
  message("  Imputing ", n_to_impute, " missing values out of ", nrow(visited), " visits")
  message("  Fixed effects: ", paste(predictors, collapse = ", "))
  message("  Random effects: intercept",
          if (!is.null(random_vars)) paste0(", ", paste(random_vars, collapse = ", ")) else "")

  # Run mice
  imp <- mice::mice(
    visited,
    m = m,
    maxit = maxit,
    predictorMatrix = pred_matrix,
    method = meth,
    where = where_mat,
    blots = blots,
    print = FALSE
  )

  # Extract completed datasets and add back unvisited rows
  datalist <- lapply(seq_len(m), function(i) {
    completed <- mice::complete(imp, action = i)
    if (nrow(unvisited) > 0) {
      rbind(completed, unvisited)
    } else {
      completed
    }
  })

  datalist
}


# ============================================
# Internal: Apply IP Method
# ============================================
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

# Null-coalescing operator (if not already defined)
`%||%` <- function(x, y) if (is.null(x)) y else x


# ============================================
# Internal: Pool Estimates
# ============================================
.pool_estimates <- function(fits, method = "mean") {

  # Extract coefficients
  coef_list <- lapply(fits, `[[`, "coef")

  # Get coefficient names (use first non-null)
  coef_names <- NULL
  for (cf in coef_list) {
    if (!is.null(names(cf))) {
      coef_names <- names(cf)
      break
    }
  }

  # Convert to matrix
  coef_mat <- do.call(rbind, coef_list)

  if (method == "mean") {
    pooled_coef <- colMeans(coef_mat, na.rm = TRUE)
    pooled_se <- apply(coef_mat, 2, sd, na.rm = TRUE) / sqrt(nrow(coef_mat))

  } else if (method == "rubin") {
    # Rubin's rules (simplified - assumes SE available)
    se_list <- lapply(fits, `[[`, "se")

    if (all(sapply(se_list, is.null))) {
      warning("SE not available; using simple averaging for variance")
      pooled_coef <- colMeans(coef_mat, na.rm = TRUE)
      pooled_se <- apply(coef_mat, 2, sd, na.rm = TRUE)
    } else {
      se_mat <- do.call(rbind, se_list)
      m <- nrow(coef_mat)

      # Within-imputation variance
      W <- colMeans(se_mat^2, na.rm = TRUE)

      # Between-imputation variance
      B <- apply(coef_mat, 2, var, na.rm = TRUE)

      # Total variance
      T_var <- W + (1 + 1/m) * B

      pooled_coef <- colMeans(coef_mat, na.rm = TRUE)
      pooled_se <- sqrt(T_var)
    }
  }

  names(pooled_coef) <- coef_names
  if (!is.null(pooled_se)) names(pooled_se) <- coef_names

  list(coef = pooled_coef, se = pooled_se)
}


# ============================================
# Print method
# ============================================
#' @export
print.MI_IP_result <- function(x, ...) {
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

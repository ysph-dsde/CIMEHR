#' Linear Increment Imputation with Informative Presence Adjustment
#'
#' @description
#' Implements a two-stage procedure that first imputes missing outcome values
#' using the linear increment (LI) model of Aalen and Gunnes (2010), then
#' applies an informative-presence (IP) adjustment method to the completed
#' dataset. The LI model is fitted with \code{\link[slim]{slim}} and
#' handles irregular observation times naturally by modelling increments
#' rather than levels.
#'
#' Unlike \code{\link{MI_IP}}, this approach produces a single completed
#' dataset (no pooling required), making it computationally faster but
#' unable to account for imputation uncertainty in standard errors.
#'
#' @param data A \code{data.frame} in long (panel) format with one row per
#'   subject-time combination. Must contain columns for subject identifier,
#'   visit time, outcome, observation indicator, and all specified
#'   covariates. Rows with \code{NA} in the time column are treated as
#'   unvisited time points.
#' @param id Character string giving the name of the subject identifier
#'   column. Default is \code{"id"}.
#' @param time Character string giving the name of the visit-time column.
#'   Default is \code{"time"}.
#' @param outcome Character string giving the name of the outcome column.
#'   Default is \code{"Y"}.
#' @param obs_indicator Character string giving the name of the binary
#'   observation indicator column (\code{1} = outcome observed, \code{0} =
#'   outcome missing at visit). Default is \code{"R"}.
#'
#' @section Linear Increment Model Parameters:
#' @param li_formula A \code{\link[stats]{formula}} for the linear increment
#'   model. The left-hand side should be the outcome variable. Default is
#'   \code{NULL}, which constructs
#'   \code{Y ~ <visit_covs> + <time>}.
#' @param li_covariance Character string specifying the covariance structure
#'   for the linear increment model. One of:
#'   \describe{
#'     \item{\code{"randomwalk"}}{(Default.) Random-walk covariance,
#'       recommended for irregular data. The variance of increments is
#'       proportional to the time gap.}
#'     \item{\code{"independence"}}{Independent increments.}
#'     \item{\code{"exchangeable"}}{Exchangeable correlation structure.}
#'   }
#' @param li_limit A one-sided \code{\link[stats]{formula}} for the limit
#'   specification passed to \code{\link[slim]{slim}}. Default is
#'   \code{~ 1} (constant limit).
#'
#' @section IP Method Parameters:
#' @param ip_method Character string specifying the IP adjustment method
#'   applied to the completed dataset. One of:
#'   \describe{
#'     \item{\code{"LMM"}}{Standard linear mixed model via
#'       \code{\link{LMM}} (no IP adjustment; included for comparison).}
#'     \item{\code{"Liang"}}{Joint model of Liang et al. (2009) via
#'       \code{\link{JMVL_IP}}.}
#'     \item{\code{"IIRR_weighting"}}{Inverse intensity rate ratio
#'       weighting via \code{\link{IIRR_balancing}} (using MLE weights).}
#'     \item{\code{"IIRR_balancing"}}{Balancing-weights estimator of
#'       Yiu and Su (2024) via \code{\link{IIRR_balancing}}.}
#'     \item{\code{"PairLik"}}{Pairwise composite likelihood via
#'       \code{\link{PairwiseLikelihood}}.}
#'   }
#' @param visit_covs Character vector of covariate column names for the
#'   visiting-process model (used by Liang and IIRR methods). Default is
#'   \code{c("Z", "X")}.
#' @param outcome_covs Character vector of covariate column names for the
#'   outcome model in the Liang method. Default is \code{NULL} (uses
#'   \code{visit_covs}).
#' @param random_covs Character vector of covariate column names with
#'   random effects in the Liang method. Default is \code{NULL} (uses the
#'   first element of \code{visit_covs}).
#' @param outcome_formula A \code{\link[stats]{formula}} for the
#'   longitudinal outcome model, used by LMM, IIRR, and PairLik methods.
#'   Default is \code{NULL}, which constructs
#'   \code{Y ~ <outcome_covs> + <time>}.
#' @param random_formula A one-sided \code{\link[stats]{formula}} for the
#'   LMM random effects (e.g., \code{~ 1 | id}). Default is \code{NULL},
#'   which constructs \code{~ 1 | <id>}.
#' @param ip_args A named list of additional arguments passed to the IP
#'   method. Common options include:
#'   \describe{
#'     \item{For IIRR methods:}{\code{phi} (numeric; sensitivity
#'       parameter, default \code{0}).}
#'     \item{For PairLik:}{\code{pair_sample} (integer; number of pairs
#'       to sample).}
#'     \item{For LMM:}{\code{visit_adjust} (\code{"OA"} or \code{"VA"}).}
#'   }
#'
#' @return An object of class \code{"LI_IP_result"}, which is a named list
#'   with the following components:
#'
#' \describe{
#'   \item{\code{coefficients}}{Named numeric vector of estimated
#'     regression coefficients from the IP method.}
#'   \item{\code{se}}{Named numeric vector of standard errors (when
#'     available from the IP method; \code{NULL} otherwise).}
#'   \item{\code{ip_method}}{Character string echoing the IP method
#'     used.}
#'   \item{\code{ip_args}}{List echoing the additional IP arguments.}
#'   \item{\code{li_model}}{A list with elements \code{formula},
#'     \code{covariance}, and \code{limit} describing the LI imputation
#'     model.}
#'   \item{\code{n_imputed}}{Integer; number of missing outcome values
#'     that were imputed.}
#'   \item{\code{imputed_data}}{A \code{data.frame}; the completed dataset
#'     after LI imputation.}
#'   \item{\code{call}}{The matched call.}
#' }
#'
#' @details
#' The procedure has two stages:
#'
#' \strong{Stage 1 --- Linear Increment Imputation.}
#' The linear increment model (Aalen and Gunnes, 2010) models the change
#' in outcome between consecutive observations:
#' \deqn{\Delta Y_{ij} = (t_j - t_{j-1})(\beta_0 + \beta^\top X_{ij})
#'   + \varepsilon_{ij},}
#' where \eqn{\text{Var}(\varepsilon_{ij}) \propto (t_j - t_{j-1})}
#' under the random-walk covariance structure. The model is fitted with
#' \code{\link[slim]{slim}} on the subset of rows with observed outcomes,
#' and predictions are generated for rows where the outcome is missing
#' but visit time is known.
#'
#' \strong{Stage 2 --- IP Adjustment.}
#' The completed dataset is passed to the selected IP method for
#' coefficient estimation. Because a single imputed dataset is used
#' (rather than multiple imputations), the reported standard errors
#' (when available) do not account for imputation uncertainty.
#'
#' S3 methods \code{print} and \code{summary} are available for objects
#' of class \code{"LI_IP_result"}.
#'
#' @references
#' Aalen, O. O. and Gunnes, N. (2010). A dynamic approach for
#' reconstructing missing longitudinal data using the linear increments
#' model. \emph{Biostatistics}, \strong{11}(3), 453--472.
#'
#' Liang, Y., Lu, W., and Ying, Z. (2009). Joint modeling and analysis of
#' longitudinal data with informative observation times. \emph{Biometrics},
#' \strong{65}(2), 377--384.
#'
#' Yiu, S. and Su, L. (2024). Accommodating informative visit times for
#' analysing irregular longitudinal data: a sensitivity analysis approach
#' with balancing weights estimators. \emph{Journal of the Royal
#' Statistical Society: Series C}, \strong{74}, 824--843.
#'
#' @seealso
#' \code{\link{MI_IP}} for the multiple imputation counterpart;
#' \code{\link[slim]{slim}} for the underlying linear increment model;
#' \code{\link{JMVL_IP}}, \code{\link{IIRR_weighting}},
#' \code{\link{IIRR_balancing}}, \code{\link{PairwiseLikelihood}}, and
#' \code{\link{LMM}} for the individual IP methods.
#'
#' @examples
#' \dontrun{
#' dat <- sim_data_gen(n = 500, seed = 42)
#'
#' # Liang's method with default LI imputation
#' result <- LI_IP(
#'   data       = dat$long_data,
#'   ip_method  = "Liang",
#'   visit_covs = c("Z", "X")
#' )
#' result
#'
#' # Custom LI model with IIRR weighting
#' result2 <- LI_IP(
#'   data          = dat$long_data,
#'   ip_method     = "IIRR_weighting",
#'   li_formula    = Y ~ Z + X + time,
#'   li_covariance = "randomwalk",
#'   visit_covs    = c("Z", "X")
#' )
#'
#' # IIRR balancing with sensitivity analysis
#' result3 <- LI_IP(
#'   data       = dat$long_data,
#'   ip_method  = "IIRR_balancing",
#'   visit_covs = c("Z", "X"),
#'   ip_args    = list(phi = 0.3)
#' )
#' summary(result3)
#' }
#'
#' @importFrom stats as.formula model.frame model.matrix predict
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
#' @noRd
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


#' Apply an Informative Presence Adjustment Method
#'
#' Internal dispatcher that routes to the appropriate IP method function.
#'
#' @param data A completed \code{data.frame} (after imputation).
#' @param ip_method Character string naming the method.
#' @param id,time,outcome,obs_indicator Column name strings.
#' @param visit_covs,outcome_covs,random_covs Character vectors of covariate
#'   names.
#' @param outcome_formula,random_formula Formula objects.
#' @param ip_args Named list of additional method-specific arguments.
#'
#' @return A list with elements \code{coef} and \code{se}.
#'
#' @keywords internal
#' @noRd
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

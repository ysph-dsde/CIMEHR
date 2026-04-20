#' Multiple Imputation (MI) with Informative Presence (IP) Adjustment
#'
#' @description
#' Implements a two-stage procedure that first imputes missing outcome values
#' using a multilevel model via \code{\link[mice]{mice}}, then applies an
#' informative-presence (IP) adjustment method to each imputed dataset and
#' pools results. The imputation model accounts for the clustered structure
#' of longitudinal data using a two-level specification (random intercept
#' and optional random slopes).
#'
#' Unlike \code{\link{LI_IP}}, this approach generates \code{m} completed
#' datasets and pools estimates, allowing imputation uncertainty to be
#' reflected in the standard errors (when Rubin's rules are used).
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
#' @param m Integer; number of multiply imputed datasets. Default is
#'   \code{5}.
#' @param maxit Integer; maximum number of iterations for
#'   \code{\link[mice]{mice}}. Default is \code{50}.
#' @param impute_method Character string specifying the \pkg{mice}
#'   imputation method for the outcome. One of:
#'   \describe{
#'     \item{\code{"2l.pan"}}{(Default.) Two-level model with heterogeneous
#'       within-group variance, implemented via \pkg{pan}.}
#'     \item{\code{"2l.lmer"}}{Two-level model using
#'       \code{\link[lme4]{lmer}}.}
#'     \item{\code{"2l.norm"}}{Two-level normal model.}
#'   }
#' @param impute_fixed Character vector of covariate column names to
#'   include as fixed effects in the imputation model. Default is
#'   \code{NULL}, which uses \code{c(visit_covs, time)}. The imputation
#'   model for \eqn{Y} is:
#'   \deqn{Y_{ij} = \beta_0 + \sum_k \beta_k X_{kij} + b_{0i} +
#'     \sum_l b_{li} Z_{lij} + \varepsilon_{ij},}
#'   where fixed effects include an intercept plus all variables in
#'   \code{impute_fixed}.
#' @param impute_random Character vector of covariate column names to
#'   include as random slopes in the imputation model. Default is
#'   \code{NULL} (random intercept only). For example, \code{c("Z")} adds
#'   a random slope for \eqn{Z}: \eqn{b_{0i} + b_{1i} Z_{ij}}.
#' @param impute_include_R Logical; whether to include the observation
#'   indicator \code{R} as a predictor in the imputation model. Default is
#'   \code{TRUE}. Recommended because \code{R} may be informative about
#'   missing \eqn{Y} values under MAR.
#'
#' @param ip_method Character string specifying the IP adjustment method
#'   applied to each completed dataset. One of:
#'   \describe{
#'     \item{\code{"LMM"}}{Standard linear mixed model via
#'       \code{\link{LMM}} (no IP adjustment; for comparison).}
#'     \item{\code{"Liang"}}{Joint model of Liang et al. (2009) via
#'       \code{\link{JMVL_Liang}}.}
#'     \item{\code{"IIRR_weighting"}}{Inverse intensity rate ratio
#'       weighting via \code{\link{IIRR_balancing}} (MLE weights).}
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
#'   LMM random effects. Default is \code{NULL}, which constructs
#'   \code{~ 1 | <id>}.
#' @param ip_args A named list of additional arguments passed to the IP
#'   method. Common options include:
#'   \describe{
#'     \item{For IIRR methods:}{\code{phi} (numeric; sensitivity parameter,
#'       default \code{0}). \code{balance_covs} (character vector).}
#'     \item{For PairLik:}{\code{pair_sample} (integer; number of pairs to
#'       sample). \code{family} (\code{"gaussian"} or
#'       \code{"binomial"}).}
#'     \item{For LMM:}{\code{method} (\code{"REML"} or \code{"ML"}).
#'       \code{optimizer} (\code{"nlminb"} or \code{"optim"}).
#'       \code{visit_adjust} (\code{NULL}, \code{"OA"}, or
#'       \code{"VA"}).}
#'   }
#'
#' @param pool_method Character string specifying the method for pooling
#'   estimates across imputations. One of:
#'   \describe{
#'     \item{\code{"mean"}}{(Default.) Simple averaging of point estimates
#'       with Monte Carlo standard errors.}
#'     \item{\code{"rubin"}}{Rubin's rules: combines within-imputation
#'       and between-imputation variance for valid inference. Requires
#'       that the IP method returns standard errors.}
#'   }
#'
#' @return An object of class \code{"MI_IP_result"}, which is a named list
#'   with the following components:
#'
#' \describe{
#'   \item{\code{coefficients}}{Named numeric vector of pooled regression
#'     coefficients.}
#'   \item{\code{se}}{Named numeric vector of pooled standard errors.}
#'   \item{\code{by_imputation}}{List of length \code{m}, each element a
#'     list with \code{coef} and (optionally) \code{se} from one imputed
#'     dataset.}
#'   \item{\code{ip_method}}{Character string echoing the IP method used.}
#'   \item{\code{ip_args}}{List echoing the additional IP arguments.}
#'   \item{\code{impute_method}}{Character string echoing the imputation
#'     method.}
#'   \item{\code{impute_model}}{A list with elements \code{fixed},
#'     \code{random}, and \code{include_R} describing the imputation
#'     model.}
#'   \item{\code{m}}{Integer; number of successful imputations.}
#'   \item{\code{m_requested}}{Integer; number of imputations originally
#'     requested.}
#'   \item{\code{call}}{The matched call.}
#' }
#'
#' @details
#' The procedure has three stages:
#'
#' \strong{Stage 1 --- Multiple Imputation.}
#' Missing outcome values (where \code{R = 0} but time is not \code{NA})
#' are imputed using a multilevel model via \code{\link[mice]{mice}}.
#' Unvisited rows (time is \code{NA}) are separated, imputation is
#' performed on the visited subset, and unvisited rows are reattached
#' afterwards. The predictor matrix is configured so that only the
#' outcome is imputed and only the specified covariates are used as
#' predictors.
#'
#' \strong{Stage 2 --- IP-Adjusted Analysis.}
#' Each completed dataset is analysed by the selected IP method. Failed
#' fits are discarded with a warning.
#'
#' \strong{Stage 3 --- Pooling.}
#' Point estimates are combined across successful imputations.
#' Under \code{pool_method = "rubin"}, Rubin's rules are applied:
#' the total variance is \eqn{T = \bar W + (1 + 1/m) B}, where
#' \eqn{\bar W} is the mean within-imputation variance and \eqn{B} is
#' the between-imputation variance.
#'
#' S3 methods \code{print} are available for objects of class
#' \code{"MI_IP_result"}.
#'
#' @references
#' Donald B. Rubin (1987). \emph{Multiple Imputation for Nonresponse in Surveys}.
#' John Wiley & Sons. \doi{10.1002/9780470316696}.
#'
#' @examples
#' dat <- sim_ehr_data[sim_ehr_data$id <= 80,
#'                     c("id", "Z", "X", "time", "R", "Y", "C")]
#' fit <- MI_IP(dat, ip_method = "LMM",
#'              visit_covs = c("Z", "X"), m = 2, maxit = 5)
#' print(fit)
#'
#' @importFrom stats as.formula sd var
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
#' @noRd
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

                     coefs <- LMM(
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
                       coef = coefs$coefficients[, "Value"],
                       se = coefs$coefficients[, "Std.Error"]
                     )
                   },

                   "Liang" = {
                     fit <- JMVL_Liang(
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

                     fit <- IIRR_balancing(
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

                     fit <- IIRR_balancing(
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

                     fit <- PairwiseLikelihood(
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

#' Null-coalescing operator
#'
#' Returns \code{y} if \code{x} is \code{NULL}, otherwise \code{x}.
#'
#' @param x,y Arbitrary objects.
#' @return \code{x} if non-\code{NULL}, else \code{y}.
#'
#' @keywords internal
#' @noRd
`%||%` <- function(x, y) if (is.null(x)) y else x


#' Pool Estimates Across Imputations
#'
#' Combines point estimates (and optionally standard errors) across multiply
#' imputed datasets using simple averaging or Rubin's rules.
#'
#' @param fits List of fit results, each with elements \code{coef} and
#'   optionally \code{se}.
#' @param method Character; \code{"mean"} or \code{"rubin"}.
#'
#' @return A list with elements \code{coef} and \code{se}.
#'
#' @keywords internal
#' @noRd
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

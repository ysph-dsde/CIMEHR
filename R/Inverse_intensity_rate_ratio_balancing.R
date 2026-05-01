#' Inverse Intensity Rate Ratio (IIRR) Estimator with Balancing Weights
#'
#' @description
#' Estimates marginal regression coefficients for a longitudinal outcome
#' using the balancing-weights approach of Yiu and Su (2024). The
#' visiting-process intensity is modelled via a partial-likelihood score
#' equation, and the outcome coefficients are obtained from weighted
#' estimating equations using two alternative weight sets: (i) balancing
#' weights that satisfy covariate-balance constraints in the spirit of
#' entropy balancing, and (ii) conventional Maximum Likelihood Estimation
#' (MLE) based inverse intensity weights.
#'
#' A sensitivity parameter \eqn{\phi} allows the analyst to explore
#' departures from the assumption that the visiting process is fully
#' captured by the observed covariates.
#'
#' @param data A \code{data.frame} in long (panel) format with one row per
#'   subject-time combination. Must contain columns for subject identifier,
#'   visit time, outcome, observation indicator, and all specified
#'   covariates.
#' @param id Character string giving the name of the subject identifier
#'   column. Default is \code{"id"}.
#' @param time Character string giving the name of the visit-time column.
#'   Default is \code{"time"}.
#' @param outcome Character string giving the name of the outcome column.
#'   Default is \code{"Y"}.
#' @param obs_indicator Character string giving the name of the binary
#'   observation indicator column (\code{1} = outcome observed at visit).
#'   Default is \code{"R"}.
#' @param visit_covs Character vector of covariate column names for the
#'   visiting-process intensity model. Default is \code{NULL}.
#' @param balance_covs Character vector of covariate column names for
#'   which covariate-balance constraints are imposed. If \code{NULL}
#'   (the default), \code{visit_covs} is used.
#' @param formula An optional \code{\link[stats]{formula}} specifying the
#'   design matrix for the longitudinal outcome model. When \code{NULL}
#'   (the default), the formula is constructed automatically as
#'   \code{outcome ~ visit_covs + time} from the column-name arguments.
#'   A user-supplied formula overrides this (e.g.,
#'   \code{Y ~ Age + Gender + time} for explicit covariate sets).
#' @param phi Numeric scalar; sensitivity parameter controlling departures
#'   from the semiparametric visiting-process model. When
#'   \eqn{\phi = 0} (the default), the estimator coincides with the
#'   standard IIRR assumption. Non-zero \eqn{\phi} tilts the balancing
#'   weights to explore sensitivity to unmeasured confounding of the
#'   visit process.
#'
#' @return An object of class \code{"Inverse_intensity_rate_ratio_balancing"}, which is a named
#'   list with the following components:
#'
#' \describe{
#'   \item{\code{beta_bal}}{Named numeric vector of outcome regression
#'     coefficients estimated with the proposed balancing weights.}
#'   \item{\code{beta_mle}}{Named numeric vector of outcome regression
#'     coefficients estimated with conventional MLE-based inverse
#'     intensity weights (for comparison).}
#'   \item{\code{gamma_hat}}{Named numeric vector of visiting-process
#'     regression coefficients, estimated by solving the partial-likelihood
#'     score equation with \code{\link[nleqslv]{nleqslv}}.}
#'   \item{\code{balancing_weights}}{Numeric vector of per-observation
#'     balancing weights.}
#'   \item{\code{mle_weights}}{Numeric vector of per-observation
#'     MLE-based inverse intensity weights.}
#'   \item{\code{breslow}}{Numeric scalar; Breslow estimate of the
#'     cumulative baseline hazard \eqn{\Lambda_0(C)}.}
#'   \item{\code{constraints}}{Numeric vector of constraint residuals
#'     (ideally near zero when balance is achieved).}
#'   \item{\code{phi}}{Numeric scalar echoing the sensitivity parameter.}
#'   \item{\code{convergence}}{A list with elements:
#'     \describe{
#'       \item{\code{gamma}}{Integer; convergence code from
#'         \code{\link[nleqslv]{nleqslv}} for the visiting-process
#'         estimation (1 = converged).}
#'       \item{\code{balancing}}{Integer; convergence code from the
#'         balancing-weights optimisation (1 = converged).}
#'     }
#'   }
#' }
#'
#' @details
#' The estimator proceeds in four steps:
#'
#' \strong{Step 1 -- Visiting-Process Coefficients.}
#' The partial-likelihood score equation for the visiting-process
#' coefficients \eqn{\gamma} is solved with
#' \code{\link[nleqslv]{nleqslv}}:
#' \deqn{\sum_j \bigl[Z_{ij} - \bar Z(\gamma; t_j)\bigr] = 0.}
#'
#' \strong{Step 2 -- MLE Weights.}
#' The MLE-based inverse intensity weight for observation \eqn{(i,j)} is
#' \deqn{w_i^{\text{MLE}} = \frac{1}{\exp(\hat\gamma^\top Z_i)
#'   \hat\Lambda_0(C)},}
#' where \eqn{\hat\Lambda_0(C)} is the Breslow cumulative baseline
#' hazard estimate.
#'
#' \strong{Step 3 -- Balancing Weights.}
#' The balancing weights \eqn{w_i^{\text{bal}}} are found by solving a
#' constrained entropy-minimisation problem:
#' \deqn{\min \sum_i w_i \log(w_i) \quad \text{s.t.} \quad
#'   \sum_i w_i X_i = \bar X,}
#' where \eqn{\bar X} is the unweighted covariate mean and the
#' constraint set is determined by \code{balance_covs}. When
#' \eqn{\phi \neq 0}, the constraints are tilted via the sensitivity
#' parameter. The dual problem is solved with
#' \code{\link[nleqslv]{nleqslv}} using Broyden's method.
#'
#' \strong{Step 4 -- Outcome Estimation.}
#' The outcome coefficients \eqn{\beta} are estimated by weighted least
#' squares using both the balancing weights (\code{beta_bal}) and the
#' MLE weights (\code{beta_mle}).
#'
#' S3 methods \code{print} and \code{summary} are available for objects
#' of class \code{"Inverse_intensity_rate_ratio_balancing"}.
#'
#' @note
#' This function requires the \pkg{nleqslv} package.
#'
#' @references
#' Sean Yiu and Li Su (2025). Accommodating informative visit times for analysing irregular
#' longitudinal data: a sensitivity analysis approach with balancing weights estimators.
#' \emph{Journal of the Royal Statistical Society Series C: Applied Statistics},
#' 74(3), 824--843. \doi{10.1093/jrsssc/qlaf002}.
#'
#' @examples
#' dat <- sim_ehr_data[sim_ehr_data$id <= 30, ]
#' fit <- Inverse_intensity_rate_ratio_balancing(dat, visit_covs = "Age")
#' print(fit)
#'
#' @importFrom stats as.formula model.frame model.matrix na.pass
#' @export
Inverse_intensity_rate_ratio_balancing <- function(data,
                                                   id = "id",
                                                   time = "time",
                                                   outcome = "Y",
                                                   obs_indicator = "R",
                                                   visit_covs = NULL,
                                                   balance_covs = NULL,
                                                   formula = NULL,
                                                   phi = 0) {

  if (!requireNamespace("nleqslv", quietly = TRUE)) {
    stop("Package 'nleqslv' is required. Install it with: install.packages('nleqslv')")
  }

  visit_covs   <- .check_covars(visit_covs,   data, "visit_covs")
  balance_covs <- .check_covars(balance_covs, data, "balance_covs")
  if (length(visit_covs) == 0L) {
    stop("`Inverse_intensity_rate_ratio_balancing` requires non-empty `visit_covs`.",
         call. = FALSE)
  }
  call <- match.call()

  # Build formula from column-name arguments if not supplied
  if (is.null(formula)) {
    formula <- as.formula(
      paste("~", paste(c(visit_covs, time), collapse = " + "))
    )
  }

  # ============================================
  # Data Preparation
  # ============================================
  data <- as.data.frame(data)
  R_vec <- data[[obs_indicator]]
  R_vec[is.na(R_vec)] <- 0
  vis_ind <- which(R_vec == 1)

  if (length(vis_ind) == 0) {
    stop("No observed visits found. Check the obs_indicator column.")
  }

  y_vec <- data[[outcome]]
  time_vec <- data[[time]]

  unique_times <- sort(unique(time_vec))
  n_obs <- nrow(data)

  # ============================================
  # Step 1: Compute selection function Q(t; phi) = exp(-phi * Y)
  # Q is only computed using Y when visited; otherwise Q = 1
  # ============================================
  Q_vec <- rep(1, n_obs)
  Q_vec[vis_ind] <- exp(-phi * y_vec[vis_ind])

  # ============================================
  # Step 2: Estimate gamma (visiting process parameters)
  # Modified score equation: sum over visits of [Z - E(Z|risk set)] * Q = 0
  # ============================================
  V_full <- as.matrix(data[, visit_covs, drop = FALSE])

  gamma_function <- function(gamma) {
    exp_Vg <- as.vector(exp(V_full %*% gamma))
    obj <- numeric(length(gamma))

    for (t_val in unique_times) {
      at_risk <- which(time_vec == t_val)
      visits_t <- intersect(vis_ind, at_risk)

      if (length(visits_t) > 0) {
        denom <- sum(exp_Vg[at_risk])
        E_Z <- colSums(V_full[at_risk, , drop = FALSE] * exp_Vg[at_risk]) / denom
        Q_visits <- Q_vec[visits_t]
        obj <- obj + colSums(V_full[visits_t, , drop = FALSE] * Q_visits) - sum(Q_visits) * E_Z
      }
    }
    return(-obj)
  }

  gamma_solve <- nleqslv::nleqslv(rep(0, length(visit_covs)), gamma_function,
                                  control = list(trace = 0, maxit = 500))
  gamma_hat <- gamma_solve$x
  names(gamma_hat) <- visit_covs

  # ============================================
  # Step 3: Compute hazard estimates = exp(gamma' * Z)
  # ============================================
  hazest <- as.vector(exp(V_full %*% gamma_hat))

  # ============================================
  # Step 4: Compute Breslow increments dLambda_0(t)
  # dLambda_0(t) = sum(Q * dN at t) / sum(hazest at risk at t)
  # ============================================
  status2 <- rep(0, n_obs)
  status2[vis_ind] <- Q_vec[vis_ind]

  breslow <- numeric(length(unique_times))
  names(breslow) <- as.character(unique_times)

  for (j in seq_along(unique_times)) {
    t_val <- unique_times[j]
    at_risk <- which(time_vec == t_val)
    numer <- sum(status2[at_risk])
    denom <- sum(hazest[at_risk])
    breslow[j] <- if (denom > 0) numer / denom else 0
  }

  # Map Breslow increments back to data
  dLambda0 <- breslow[match(as.character(time_vec), names(breslow))]
  dLambda0[is.na(dLambda0)] <- 0

  # ============================================
  # Step 5: Build design matrix for balancing
  # Includes: intercept, time, covariates, covariate*time interactions
  # ============================================
  # `.check_covars()` returns character(0) for NULL, so use a length check
  # rather than is.null(). By default, balance on the same covariates used
  # in the visiting-process model.
  if (length(balance_covs) == 0L) balance_covs <- visit_covs

  V_bal <- as.matrix(data[, balance_covs, drop = FALSE])

  DesignMat <- cbind(
    intercept = 1,
    time = time_vec,
    V_bal,
    V_bal * time_vec
  )

  colnames(DesignMat) <- c(
    "intercept",
    "time",
    balance_covs,
    paste0(balance_covs, "_x_time")
  )

  # ============================================
  # Step 6: Compute balancing constraints (target values)
  # Target: sum of h(t) weighted by dLambda_0(t) across entire population
  # ============================================
  constrain <- colSums(DesignMat * dLambda0)

  # ============================================
  # Step 7: Estimate balancing weights
  # Solve: sum(h * W) = constrain, where W = exp(gamma_b' * h) * Q
  # ============================================
  DesignMat_vis <- DesignMat[vis_ind, , drop = FALSE]
  offset_vis <- Q_vec[vis_ind]

  bal_fit_fun_sa <- function(MAT, cons, offs) {
    bal_fun <- function(gamma) {
      wei <- as.vector(exp(MAT %*% gamma)) * offs
      colSums(MAT * wei) - cons
    }

    jac_fun_bal <- function(gamma) {
      wei <- as.vector(exp(MAT %*% gamma)) * offs
      t(MAT) %*% (MAT * wei)
    }

    result <- nleqslv::nleqslv(rep(0, ncol(MAT)), bal_fun, jac = jac_fun_bal,
                               method = "Broyden", control = list(maxit = 1000))

    list(
      weights = as.vector(exp(MAT %*% result$x)) * offs,
      convergence = result$termcd
    )
  }

  bal_result <- tryCatch({
    bal_fit_fun_sa(DesignMat_vis, constrain, offset_vis)
  }, error = function(e) {
    msg <- conditionMessage(e)
    msg <- sub("\\.*$", "", msg)
    message(
      "Note: ", msg,
      ". Check initial x and/or correctness of function. ",
      "Falling back to Maximum Likelihood Estimation (MLE) weights."
    )
    list(
      weights = (1 / hazest[vis_ind]) * offset_vis,
      convergence = -1
    )
  })

  Bal_weights <- bal_result$weights
  bal_convergence <- bal_result$convergence

  # ============================================
  # Step 8: Compute MLE weights (for comparison)
  # MLE weights = (1/hazest) * Q
  # ============================================
  MLE_weights <- (1 / hazest[vis_ind]) * offset_vis

  # ============================================
  # Step 9: Estimate beta using weights
  # Weighted least squares: (X'WX)^{-1} X'WY
  # ============================================
  data_obs <- data[vis_ind, ]
  y_obs <- y_vec[vis_ind]

  mf <- model.frame(formula, data = data_obs, na.action = na.pass)
  X <- model.matrix(formula, data = mf)

  # Balancing weights estimate

  beta_bal <- tryCatch({
    as.vector(solve(t(X) %*% (X * Bal_weights)) %*% (t(X) %*% (Bal_weights * y_obs)))
  }, error = function(e) {
    warning("Beta estimation with balancing weights failed: ", e$message)
    rep(NA, ncol(X))
  })
  names(beta_bal) <- colnames(X)

  # MLE weights estimate
  beta_mle <- tryCatch({
    as.vector(solve(t(X) %*% (X * MLE_weights)) %*% (t(X) %*% (MLE_weights * y_obs)))
  }, error = function(e) {
    warning("Beta estimation with MLE weights failed: ", e$message)
    rep(NA, ncol(X))
  })
  names(beta_mle) <- colnames(X)

  # ============================================
  # Return results
  # ============================================
  return(structure(list(
    beta_bal = beta_bal,
    beta_mle = beta_mle,
    gamma_hat = gamma_hat,
    balancing_weights = Bal_weights,
    mle_weights = MLE_weights,
    breslow = breslow,
    constraints = constrain,
    phi = phi,
    convergence = list(
      gamma = gamma_solve$termcd,
      balancing = bal_convergence
    ),
    call = call
  ), class = "Inverse_intensity_rate_ratio_balancing"))
}

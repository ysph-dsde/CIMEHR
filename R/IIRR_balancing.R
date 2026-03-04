#' Balancing Weights Inverse Intensity Rate Ratio Estimator
#'
#' Implements the balancing weights approach from Yiu and Su (2024) for
#' accommodating informative visit times in longitudinal data analysis.
#' This function estimates marginal regression coefficients using both

#' balancing weights (proposed method) and MLE weights (standard method)
#' for comparison.
#'
#' @param data Data frame in panel format (id, time, Y, R, covariates).
#' @param id Character. Subject ID column name. Default is "id".
#' @param time Character. Time variable column name. Default is "time".
#' @param outcome Character. Outcome variable column name. Default is "Y".
#' @param obs_indicator Character. Observation indicator column name
#'        (1 = visited/observed, 0 = not visited). Default is "R".
#' @param visit_covs Character vector. Covariate names for visiting process model.
#'        Default is c("Z", "X").
#' @param balance_covs Character vector. Covariates to balance. Default is NULL,
#'        which uses the same covariates as visit_covs. Time interactions are
#'        automatically added.
#' @param formula Formula. Model formula for longitudinal outcome model.
#'        Default is ~ Z + X + time.
#' @param phi Numeric. Sensitivity parameter for selection function. Default is 0.
#'        - phi = 0: Visiting at Random (VAR) assumption
#'        - phi > 0: Higher outcome increases visit probability
#'        - phi < 0: Higher outcome decreases visit probability
#'
#' @return A list containing:
#'   \item{beta_bal}{Numeric vector. Outcome model coefficients estimated using
#'         balancing weights. This is the proposed method from Yiu and Su (2024),
#'         which provides more robust and efficient estimates than MLE weights.}
#'   \item{beta_mle}{Numeric vector. Outcome model coefficients estimated using
#'         standard MLE-based inverse intensity weights. Provided for comparison
#'         with the balancing weights approach.}
#'   \item{gamma_hat}{Numeric vector. Visiting process coefficients from the
#'         Cox-type model for visit intensity.}
#'   \item{balancing_weights}{Numeric vector. Estimated balancing weights for
#'         each observed visit. These weights satisfy the balancing constraints
#'         ensuring covariate distributions among visitors match the target population.
#'         Length equals the number of visits (sum of obs_indicator == 1).}
#'   \item{mle_weights}{Numeric vector. Standard inverse intensity weights
#'         calculated as (1/hazest) * Q. Length equals the number of visits.}
#'   \item{breslow}{Numeric vector. Breslow increments of baseline hazard dLambda_0(t)
#'         at each unique time point.}
#'   \item{constraints}{Numeric vector. Target values for balancing constraints.
#'         Represents covariate means in the target population weighted by baseline hazard.}
#'   \item{phi}{Numeric. The sensitivity parameter used in the analysis.}
#'   \item{convergence}{List containing convergence codes for gamma estimation
#'         and balancing weight estimation. Code 1 indicates successful convergence.}
#'
#' @details
#' The function implements two weighting approaches:
#'
#' **Balancing Weights (beta_bal)**: Weights are estimated by solving balancing
#' equations that ensure the weighted covariate distribution among visitors matches
#' the distribution in the target population. This approach is more robust to model
#' misspecification and has smaller finite-sample bias compared to MLE weights.
#'
#' **MLE Weights (beta_mle)**: Standard inverse intensity weights estimated from
#' the Cox partial likelihood. May have larger bias when the visit process is
#' highly dependent on covariates.
#'
#' Both methods can incorporate the sensitivity parameter phi to assess the impact
#' of informative visit times on the results.
#'
#' @references
#' Yiu, S. and Su, L. (2024). Accommodating informative visit times for analysing
#' irregular longitudinal data: a sensitivity analysis approach with balancing
#' weights estimators. Journal of the Royal Statistical Society Series C:
#' Applied Statistics, 74, 824-843.
#'
#' @examples
#' \dontrun{
#' # Basic usage with VAR assumption (phi = 0)
#' result <- IIRR_balancing(data, visit_covs = c("Z", "X"), phi = 0)
#'
#' # Sensitivity analysis with phi = 0.3
#' result_sa <- IIRR_balancing(data, visit_covs = c("Z", "X"), phi = 0.3)
#'
#' # Compare balancing vs MLE weights
#' result$beta_bal  # Balancing weights (recommended)
#' result$beta_mle  # MLE weights (for comparison)
#' }
#'
#' @export
IIRR_balancing <- function(data,
                               id = "id",
                               time = "time",
                               outcome = "Y",
                               obs_indicator = "R",
                               visit_covs = c("Z", "X"),
                               balance_covs = NULL,
                               formula = ~ Z + X + time,
                               phi = 0) {

  if (!requireNamespace("nleqslv", quietly = TRUE)) {
    stop("Package 'nleqslv' is required. Install it with: install.packages('nleqslv')")
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
  id_vec <- data[[id]]

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
  if (is.null(balance_covs)) balance_covs <- visit_covs

  DesignMat <- cbind(
    intercept = 1,
    time = time_vec,
    V_full,
    V_full * time_vec
  )
  colnames(DesignMat) <- c("intercept", "time", balance_covs,
                           paste0(balance_covs, "_x_time"))

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
    warning("Balancing weight estimation failed: ", e$message,
            ". Falling back to MLE weights.")
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
  return(list(
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
    )
  ))
}

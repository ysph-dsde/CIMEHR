#' Inverse Intensity Rate Ratio Weighting Estimator
#'
#' Fits a weighted GEE model using inverse intensity weights from the visiting process.
#' Only uses observed data (R == 1 or R is NA).
#'
#' @param data Data frame containing longitudinal data.
#' @param id Subject ID column name. Default is "id".
#' @param time Time variable column name. Default is "time".
#' @param outcome Outcome variable column name. Default is "Y".
#' @param obs_indicator Observation indicator column name. Default is "R".
#' @param visit_covs Character vector of covariate names for visiting process. Default is c("Z", "X").
#' @param formula Formula for longitudinal model (e.g., ~ Z + X + time). Default is ~ Z + X + time.
#'
#' @return A list containing:
#'   \item{gamma_hat}{Visiting process coefficients}
#'   \item{beta_hat}{Longitudinal process coefficients}
#'
#' @export
IIRR_weighting <- function(data,
                               id = "id",
                               time = "time",
                               outcome = "Y",
                               obs_indicator = "R",
                               visit_covs = c("Z", "X"),
                               formula = ~ Z + X + time) {

  if (!requireNamespace("nleqslv", quietly = TRUE)) install.packages("nleqslv")

  # Filter to observed data only (R == 1 or R is NA)
  R_vec <- data[[obs_indicator]]
  data <- data[R_vec == 1 | is.na(R_vec), ]

  # Visiting process variables
  V <- as.matrix(data[, c(id, visit_covs)])
  VV <- as.matrix(V[!duplicated(V), -1])
  V <- as.matrix(data[, visit_covs])

  y <- data[[outcome]]

  # ============================================
  # Step 1: Estimate gamma (visiting process)
  # ============================================
  gamma_function <- function(gamma) {
    gamma <- matrix(gamma, ncol = 1)
    V_denom <- sum(exp(VV %*% gamma))
    V_numer <- colSums(VV * as.numeric(exp(VV %*% gamma)))
    obj <- colSums(V) - (V_numer / V_denom) * nrow(V)
    return(obj)
  }

  gamma_init <- rep(1, ncol(V))
  gamma_solve <- nleqslv::nleqslv(x = gamma_init, fn = gamma_function, control = list(trace = 0))
  gamma_hat <- gamma_solve$x
  names(gamma_hat) <- visit_covs

  # ============================================
  # Step 2: Compute weights
  # ============================================
  weights <- as.vector(exp(V %*% gamma_hat))

  # ============================================
  # Step 3: Build design matrix from formula
  # ============================================
  mf <- model.frame(formula, data = data, na.action = na.pass)
  X_GEE <- model.matrix(formula, data = mf)

  # ============================================
  # Step 4: Estimate beta (weighted estimating equation)
  # ============================================
  Y_function <- function(beta) {
    res <- as.vector((1 / weights) * (y - X_GEE %*% beta))
    obj <- colSums(na.omit(X_GEE * res))
    return(obj)
  }

  beta_init <- rep(0.1, ncol(X_GEE))
  beta_solve <- nleqslv::nleqslv(x = beta_init, fn = Y_function, control = list(trace = 0))
  beta_hat <- beta_solve$x
  names(beta_hat) <- colnames(X_GEE)

  return(list(
    gamma_hat = gamma_hat,
    beta_hat = beta_hat
  ))
}

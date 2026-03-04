#' EHR Joint Model for Informative Visiting Process
#'
#' Fits a joint model accounting for informative visiting and observation processes in EHR data.
#'
#' @param data Data frame containing longitudinal data.
#' @param id Subject ID column name. Default is "id".
#' @param time Time variable column name. Default is "time".
#' @param outcome Outcome variable column name. Default is "Y".
#' @param obs_indicator Observation indicator column name (1 = observed, 0 = missing). Default is "R".
#' @param visit_covs Character vector of covariate names for visiting process. Default is c("Z", "X").
#' @param long_covs Character vector of covariate names for longitudinal process. Default is c("Z", "X").
#' @param random_covs Character vector of covariate names with random effects. Default is "Z".
#'
#' @return A list containing:
#'   \item{gamma_hat}{Visiting process coefficients}
#'   \item{alpha_hat}{Observation process coefficients}
#'   \item{Lambda_C_est}{Estimated cumulative intensity}
#'   \item{sigma_sq_hat}{Estimated random effect variance}
#'   \item{beta_hat}{Longitudinal process coefficients}
#'   \item{theta_hat}{Random effect loading coefficients}
EHRJoint <- function(data,
                              id = "id",
                              time = "time",
                              outcome = "Y",
                              obs_indicator = "R",
                              visit_covs = c("Z", "X"),
                              long_covs = c("Z", "X"),
                              random_covs = "Z") {

  if (!requireNamespace("nleqslv", quietly = TRUE)) install.packages("nleqslv")

  # Visiting process variables
  V <- as.matrix(data[, c(id, visit_covs)])
  VV <- as.matrix(V[!duplicated(V), -1])
  V <- as.matrix(data[, visit_covs])

  # Longitudinal process variables
  X <- as.matrix(data[, c(id, long_covs)])
  XX <- as.matrix(X[!duplicated(X), -1])
  X <- as.matrix(data[, long_covs])

  # Random effect variables
  Q <- data[, random_covs, drop = FALSE]

  # Compute ni
  id_vec <- data[[id]]
  unique_ids <- unique(id_vec)
  R <- data[[obs_indicator]]

  ni <- sapply(unique_ids, function(i) sum(!is.na(R[id_vec == i])))
  ni_0to1 <- ifelse(ni == 0, 1, ni)

  Nit <- as.numeric(!is.na(data[[time]]))
  C <- rep(max(data[[time]], na.rm = TRUE), nrow(VV))
  Time <- data[[time]]
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
  # Step 2: Estimate alpha (observation process)
  # ============================================
  obs_formula <- as.formula(paste(obs_indicator, "~", paste(visit_covs, collapse = " + ")))
  fit_IP <- glm(obs_formula, data = data, family = binomial)
  alpha_hat <- coef(fit_IP)

  wi <- as.vector(plogis(cbind(1, as.matrix(data[, visit_covs])) %*% alpha_hat))
  wi_base <- wi[!duplicated(id_vec)]

  # ============================================
  # Step 3: Estimate Lambda_C
  # ============================================
  exp_gamma <- as.vector(exp(VV %*% gamma_hat))
  Lambda_C_est <- sum(!is.na(Time)) / sum(exp_gamma)

  # ============================================
  # Step 4: Estimate sigma^2
  # ============================================
  exp_gamma_Lambda <- exp_gamma * Lambda_C_est
  sigma_sq_hat <- max((sum(ni^2 - ni - exp_gamma_Lambda^2) / sum(exp_gamma_Lambda^2)), 0)

  # ============================================
  # Step 5: Estimate B
  # ============================================
  Q0 <- cbind(1, Q)
  B_i <- ((1 + ni * sigma_sq_hat) / (1 + exp_gamma_Lambda * sigma_sq_hat) - 1)
  Bhat_long <- rep(B_i, ni_0to1)
  Bhat <- Bhat_long * Q0
  Bhat_base0 <- cbind(id_vec, Bhat)
  Bhat_base <- Bhat_base0[!duplicated(id_vec), -1]

  # ============================================
  # Step 6: Estimate beta
  # ============================================
  denom <- rep(sum(ni * wi_base / Lambda_C_est), length(Time))

  numer_X <- apply(XX, 2, function(x) {
    rep(sum(ni * x * wi_base / Lambda_C_est), length(Time))
  })
  numer_B <- apply(Bhat_base, 2, function(x) {
    rep(sum(ni * x * wi_base / Lambda_C_est, na.rm = TRUE), length(Time))
  })

  if (sigma_sq_hat != 0) {
    X_bar <- numer_X / as.vector(denom)
    B_bar <- numer_B / as.vector(denom)
    design_M <- as.matrix(cbind(X, Bhat))
    design_Mbar <- as.matrix(cbind(X_bar, B_bar))

    Liang_function <- function(beta) {
      U_Liang <- (design_M - design_Mbar) * as.vector(y - design_M %*% beta)
      obj <- colSums(na.omit(U_Liang))
      return(obj)
    }

    beta_init <- rep(0.1, ncol(design_M))
    beta_solve <- nleqslv::nleqslv(x = beta_init, fn = Liang_function, control = list(trace = 0))
    beta_hat <- beta_solve$x[1:ncol(X)]
    theta_hat <- beta_solve$x[-(1:ncol(X))]

  } else {
    X_bar <- numer_X / denom
    design_M <- as.matrix(X)
    design_Mbar <- as.matrix(numer_X / denom)

    Liang_function <- function(beta) {
      U_Liang <- (design_M - design_Mbar) * as.vector(y - design_M %*% beta)
      obj <- colSums(na.omit(U_Liang))
      return(obj)
    }

    beta_init <- rep(0, ncol(design_M))
    beta_solve <- nleqslv::nleqslv(x = beta_init, fn = Liang_function, control = list(trace = 0))
    beta_hat <- beta_solve$x[1:ncol(V)]
    theta_hat <- c(NA, NA)
  }

  names(beta_hat) <- long_covs

  return(list(
    gamma_hat = gamma_hat,
    alpha_hat = alpha_hat,
    Lambda_C_est = Lambda_C_est,
    sigma_sq_hat = sigma_sq_hat,
    beta_hat = beta_hat,
    theta_hat = theta_hat
  ))
}

#fit2 = EHRJoint(dat$long_data)

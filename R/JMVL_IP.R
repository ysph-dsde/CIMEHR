#' Joint Model for Visiting and Longitudinal Processes (Liang's Method)
#'
#' @description
#' Fits a joint model for longitudinal data with an informative visiting
#' process using the estimating-equations approach of Liang, Lu, and Ying
#' (2009). The visiting-process intensity is modelled via a Cox-type
#' partial-likelihood score equation, and the longitudinal outcome is
#' estimated by solving Liang-type risk-set-centred estimating equations
#' that incorporate an empirical Bayes estimate of the shared latent
#' factor.
#'
#' Unlike \code{\link{EHRJoint}}, this method does not model a separate
#' observation process; it assumes that the outcome is always observed at
#' visits (or equivalently, filters to rows with \code{R == 1} or
#' \code{R} is \code{NA}).
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
#'   observation indicator column. Rows with \code{R == 1} or \code{R} is
#'   \code{NA} are retained; all others are dropped. Default is
#'   \code{"R"}.
#' @param visit_covs Character vector of covariate column names for the
#'   visiting-process intensity model. Default is \code{c("Z", "X")}.
#' @param long_covs Character vector of covariate column names for the
#'   fixed-effect part of the longitudinal outcome model. Default is
#'   \code{c("Z", "X")}.
#' @param random_covs Character vector of covariate column names that
#'   interact with the latent factor in the outcome model. An intercept
#'   column is prepended automatically. Default is \code{"Z"}.
#'
#' @return A named list with four components:
#'
#' \describe{
#'   \item{\code{gamma_hat}}{Named numeric vector of visiting-process
#'     regression coefficients, estimated by solving the partial-likelihood
#'     score equation with \code{\link[nleqslv]{nleqslv}}.}
#'   \item{\code{Lambda_C_est}}{Numeric scalar; Breslow estimate of the
#'     cumulative baseline hazard \eqn{\Lambda_0(C)} over the full
#'     follow-up window.}
#'   \item{\code{sigma_sq_hat}}{Numeric scalar; method-of-moments estimate
#'     of the frailty variance \eqn{\sigma^2}. Constrained to be
#'     non-negative.}
#'   \item{\code{beta_hat}}{Named numeric vector of outcome fixed-effect
#'     coefficients, estimated from the Liang-type estimating equations.}
#' }
#'
#' @details
#' The model consists of two linked components:
#'
#' \strong{Visiting Process.}
#' The visit intensity is
#' \deqn{\lambda_i(t) = \eta_i \exp(\gamma^\top Z_i)\,\lambda_0(t),}
#' where \eqn{\eta_i} is a shared frailty. The coefficient vector
#' \eqn{\gamma} is estimated by solving the partial-likelihood score
#' equation, the cumulative baseline hazard \eqn{\Lambda_0(C)} by the
#' Breslow estimator, and the frailty variance \eqn{\sigma^2} by method
#' of moments.
#'
#' \strong{Outcome Model.}
#' An empirical Bayes estimator of the latent loading \eqn{B_i} is
#' constructed as
#' \deqn{\hat B_i = \frac{1 + n_i \hat\sigma^2}{1 + \hat\Lambda_i
#'   \hat\sigma^2} - 1,}
#' and the outcome model
#' \deqn{Y_i(t) = \beta^\top X_i(t) + B_i Z_i^*(t) + \varepsilon_i(t)}
#' is estimated by solving risk-set-centred estimating equations
#' \deqn{\sum_{i,j} \bigl[(M_{ij} - \bar M(t_j))\bigr]
#'   (Y_{ij} - M_{ij}^\top \beta^*) = 0,}
#' where \eqn{M_{ij} = (X_{ij}, \hat B_i Z_{ij}^*)} and
#' \eqn{\bar M(t)} is the risk-set-weighted mean.
#'
#' This differs from \code{\link{EHRJoint}} in that no observation-process
#' weights enter the estimating equations.
#'
#' @note
#' This function requires the \pkg{nleqslv} package. It will attempt to
#' install the package automatically if it is not available.
#'
#' No standard errors are returned. For inference, consider bootstrapping
#' the entire procedure.
#'
#' @references
#' Yu Liang, Wenbin Lu, and Zhiliang Ying (2009). Joint modeling and analysis of longitudinal
#' data with informative observation times. \emph{Biometrics}, 65(2), 377--384.
#' \doi{10.1111/j.1541-0420.2008.01104.x}.
#'
#' @examples
#' \dontrun{
#' dat <- sim_data_gen(n = 500, seed = 42)
#'
#' fit <- JMVL_IP(
#'   data       = dat$long_data,
#'   visit_covs = c("Z", "X"),
#'   long_covs  = c("Z", "X"),
#'   random_covs = "Z"
#' )
#' fit$gamma_hat
#' fit$beta_hat
#' fit$sigma_sq_hat
#' }
#'
#' @export
JMVL_IP <- function(data,
                    id = "id",
                    time = "time",
                    outcome = "Y",
                    obs_indicator = "R",
                    visit_covs = c("Z", "X"),
                    long_covs = c("Z", "X"),
                    random_covs = "Z") {

  if (!requireNamespace("nleqslv", quietly = TRUE)) {
    stop("Package 'nleqslv' is required. Install it with: install.packages('nleqslv')",
         call. = FALSE)
  }

  # Filter to observed data only (R == 1 or R is NA)
  R_vec <- data[[obs_indicator]]
  data <- data[R_vec == 1 | is.na(R_vec), ]

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
  time_vec <- data[[time]]

  ni <- sapply(unique_ids, function(i) sum(!is.na(time_vec[id_vec == i])))
  ni_0to1 <- ifelse(ni == 0, 1, ni)

  C <- rep(max(time_vec, na.rm = TRUE), nrow(VV))
  Time <- time_vec
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
  # Step 2: Estimate Lambda_C
  # ============================================
  exp_gamma <- as.vector(exp(VV %*% gamma_hat))
  Lambda_C_est <- sum(!is.na(Time)) / sum(exp_gamma)

  # ============================================
  # Step 3: Estimate sigma^2
  # ============================================
  exp_gamma_Lambda <- exp_gamma * Lambda_C_est
  sigma_sq_hat <- max((sum(ni^2 - ni - exp_gamma_Lambda^2) / sum(exp_gamma_Lambda^2)), 0)

  # ============================================
  # Step 4: Estimate B
  # ============================================
  Q0 <- cbind(1, Q)
  B_i <- ((1 + ni * sigma_sq_hat) / (1 + exp_gamma_Lambda * sigma_sq_hat) - 1)
  Bhat_long <- rep(B_i, ni_0to1)
  Bhat <- Bhat_long * Q0
  Bhat_base0 <- cbind(id_vec, Bhat)
  Bhat_base <- Bhat_base0[!duplicated(id_vec), -1]

  # ============================================
  # Step 5: Estimate beta (no wi weighting)
  # ============================================
  denom <- rep(sum(ni / Lambda_C_est), length(Time))

  numer_X <- apply(XX, 2, function(x) {
    rep(sum(ni * x / Lambda_C_est), length(Time))
  })
  numer_B <- apply(Bhat_base, 2, function(x) {
    rep(sum(ni * x / Lambda_C_est, na.rm = TRUE), length(Time))
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
  }

  names(beta_hat) <- long_covs

  return(list(
    gamma_hat = gamma_hat,
    Lambda_C_est = Lambda_C_est,
    sigma_sq_hat = sigma_sq_hat,
    beta_hat = beta_hat
  ))
}

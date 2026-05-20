#' Joint Model for Visiting and Longitudinal (JMVL) Processes with
#' Informative Presence (IP) Adjustment (Liang's Method)
#'
#' @description
#' Fits a joint model for longitudinal data with an informative visiting
#' process using the estimating-equations approach of Liang, Lu, and Ying
#' (2009). The visiting-process intensity is modelled via a Cox-type
#' partial-likelihood score equation, and the longitudinal outcome is
#' estimated by solving Liang-type risk-set-centred estimating equations
#' that incorporate an empirical Bayes estimate of the shared latent
#' factor. Analytic standard errors are returned for both the visiting
#' process (model-based, from the inverse information matrix) and the
#' outcome (clustered sandwich, conditional on the Stage 1 nuisance
#' estimates).
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
#' @return An object of class
#'   \code{"Joint_modeling_visiting_and_longitudinal_Liang"}, a named list
#'   with the following components:
#'
#' \describe{
#'   \item{\code{gamma_hat}}{Named numeric vector of visiting-process
#'     regression coefficients, estimated by solving the partial-likelihood
#'     score equation with \code{\link[nleqslv]{nleqslv}}.}
#'   \item{\code{se_gamma}}{Named numeric vector of model-based standard
#'     errors for \code{gamma_hat}, obtained from the inverse observed
#'     information matrix.}
#'   \item{\code{vcov_gamma}}{Numeric matrix; the model-based
#'     variance-covariance matrix for \code{gamma_hat}.}
#'   \item{\code{Lambda_C_est}}{Numeric scalar; Breslow estimate of the
#'     cumulative baseline hazard \eqn{\Lambda_0(C)} over the full
#'     follow-up window.}
#'   \item{\code{sigma_sq_hat}}{Numeric scalar; method-of-moments estimate
#'     of the frailty variance \eqn{\sigma^2}. Constrained to be
#'     non-negative.}
#'   \item{\code{beta_hat}}{Named numeric vector of outcome fixed-effect
#'     coefficients, estimated from the Liang-type estimating equations.}
#'   \item{\code{se_beta}}{Named numeric vector of clustered sandwich
#'     standard errors for \code{beta_hat}. The sandwich estimator is
#'     conditional on the Stage 1 estimates
#'     \eqn{(\hat\gamma, \hat\sigma^2, \hat\Lambda_0)}.}
#'   \item{\code{theta_hat}}{Named numeric vector of latent-link
#'     coefficients (the part of the augmented coefficient vector
#'     multiplying \eqn{\hat B_i Z^*_{ij}}). \code{NULL} when
#'     \eqn{\hat\sigma^2 = 0}.}
#'   \item{\code{se_theta}}{Named numeric vector of clustered sandwich
#'     standard errors for \code{theta_hat}. \code{NULL} when
#'     \eqn{\hat\sigma^2 = 0}.}
#'   \item{\code{vcov_outcome}}{Numeric matrix; the full clustered
#'     sandwich variance-covariance matrix for the augmented coefficient
#'     vector \code{c(beta_hat, theta_hat)}.}
#'   \item{\code{call}}{The matched call.}
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
#' of moments. Standard errors come from the inverse observed information
#' matrix.
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
#' \eqn{\bar M(t)} is the risk-set-weighted mean. Because this
#' equation is linear in \eqn{\beta^*}, it admits a closed-form
#' solution \eqn{\hat\beta^* = (H^\top M)^{-1} H^\top Y} with
#' \eqn{H = M - \bar M}.
#'
#' \strong{Sandwich Variance.}
#' The clustered sandwich estimator for \eqn{\hat\beta^*} is
#' \deqn{\hat V = B^{-1}\Bigl(\sum_i A_i A_i^\top\Bigr) B^{-\top},
#'   \qquad B = H^\top M,\quad
#'   A_i = \sum_j H_{ij}(Y_{ij} - M_{ij}^\top\hat\beta^*).}
#' This variance is conditional on the Stage 1 nuisance estimates;
#' it does not propagate estimation uncertainty from
#' \eqn{(\hat\gamma, \hat\sigma^2, \hat\Lambda_0)}.
#'
#' This differs from \code{\link{EHRJoint}} in that no observation-process
#' weights enter the estimating equations.
#'
#' @note
#' This function requires the \pkg{nleqslv} and \pkg{MASS} packages.
#'
#' @references
#' Liang, Y., Lu, W., and Ying, Z. (2009). Joint modeling and analysis of
#' longitudinal data with informative observation times.
#' \emph{Biometrics}, 65(2), 377--384.
#' \doi{10.1111/j.1541-0420.2008.01104.x}.
#'
#' @seealso
#' \code{\link{EHRJoint}} for a three-stage joint model that additionally
#' incorporates an observation-process component;
#' \code{\link{CIMEHR}} for the three-stage CIMEHR estimator with a full
#' shared latent structure.
#'
#' @examples
#' dat <- sim_ehr_data[sim_ehr_data$id <= 30, ]
#' fit <- Joint_modeling_visiting_and_longitudinal_Liang(
#'   dat,
#'   visit_covs  = "Age",
#'   long_covs   = "Age",
#'   random_covs = "Age"
#' )
#' print(fit)
#'
#' @importFrom stats complete.cases
#' @export
Joint_modeling_visiting_and_longitudinal_Liang <- function(data,
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
  if (!requireNamespace("MASS", quietly = TRUE)) {
    stop("Package 'MASS' is required. Install it with: install.packages('MASS')",
         call. = FALSE)
  }

  call <- match.call()

  visit_covs  <- .check_covars(visit_covs,  data, "visit_covs")
  long_covs   <- .check_covars(long_covs,   data, "long_covs")
  random_covs <- .check_covars(random_covs, data, "random_covs")
  if (length(visit_covs)  == 0L) stop("`visit_covs` must include at least one covariate.",  call. = FALSE)
  if (length(long_covs)   == 0L) stop("`long_covs` must include at least one covariate.",   call. = FALSE)
  if (length(random_covs) == 0L) stop("`random_covs` must include at least one covariate.", call. = FALSE)

  # ============================================
  # Data prep
  # ============================================

  # Filter to observed data only (R == 1 or R is NA)
  R_vec <- data[[obs_indicator]]
  data  <- data[R_vec == 1 | is.na(R_vec), ]
  data  <- data[order(data[[id]], data[[time]]), ]

  id_vec     <- data[[id]]
  unique_ids <- unique(id_vec)
  n_subj     <- length(unique_ids)
  time_vec   <- data[[time]]
  y          <- data[[outcome]]

  sf <- which(!duplicated(id_vec))    # first-row indices per subject

  # Visiting-process design matrices
  V_full <- as.matrix(data[, visit_covs, drop = FALSE])   # all rows
  VV     <- V_full[sf, , drop = FALSE]                    # one row per subject

  # Outcome design matrices
  X  <- as.matrix(data[, long_covs, drop = FALSE])        # all rows
  XX <- X[sf, , drop = FALSE]                             # one row per subject

  # Random-effect design (intercept prepended)
  Q  <- as.matrix(data[, random_covs, drop = FALSE])
  Q0 <- cbind(Intercept = 1, Q)

  # Per-subject visit count
  ni      <- tabulate(match(id_vec, unique_ids), nbins = n_subj)
  ni_0to1 <- pmax(ni, 1L)

  # ============================================
  # Step 1: Estimate gamma (visiting process)
  # ============================================

  gamma_function <- function(gamma) {
    lin_pred <- as.vector(VV %*% gamma)
    w     <- exp(lin_pred - max(lin_pred))      # log-sum-exp for stability
    w_sum <- sum(w)
    colSums(V_full) - (colSums(VV * w) / w_sum) * nrow(V_full)
  }

  gamma_init  <- rep(0, ncol(V_full))
  gamma_solve <- nleqslv::nleqslv(x = gamma_init, fn = gamma_function,
                                  control = list(trace = 0))
  gamma_hat <- gamma_solve$x
  names(gamma_hat) <- visit_covs

  # Model-based SE for gamma (inverse information)
  exp_gW         <- as.vector(exp(VV %*% gamma_hat))
  sum_risk       <- sum(exp_gW)
  Xbar_global    <- colSums(exp_gW * VV) / sum_risk
  VV_centered    <- sweep(VV, 2, Xbar_global, "-")
  Info_gamma     <- nrow(V_full) * crossprod(VV_centered * sqrt(exp_gW)) / sum_risk
  inv_Info_gamma <- tryCatch(solve(Info_gamma),
                             error = function(e) MASS::ginv(Info_gamma))
  vcov_gamma     <- inv_Info_gamma
  se_gamma       <- sqrt(pmax(diag(vcov_gamma), 0))
  names(se_gamma) <- visit_covs
  dimnames(vcov_gamma) <- list(visit_covs, visit_covs)

  # ============================================
  # Step 2: Estimate Lambda_C (Breslow)
  # ============================================

  exp_gamma    <- as.vector(exp(VV %*% gamma_hat))
  Lambda_C_est <- sum(!is.na(time_vec)) / sum(exp_gamma)

  # ============================================
  # Step 3: Estimate sigma^2 (method of moments)
  # ============================================

  exp_gamma_Lambda <- exp_gamma * Lambda_C_est
  sigma_sq_hat     <- max(sum(ni^2 - ni - exp_gamma_Lambda^2) /
                            sum(exp_gamma_Lambda^2), 0)

  # ============================================
  # Step 4: Construct B_hat (empirical Bayes)
  # ============================================

  B_i  <- (1 + ni * sigma_sq_hat) / (1 + exp_gamma_Lambda * sigma_sq_hat) - 1
  Bhat <- rep(B_i, ni_0to1) * Q0

  Bhat_subj <- Bhat[sf, , drop = FALSE]   # one row per subject (for M_bar)

  # ============================================
  # Step 5: Estimate beta (direct linear solve)
  # ============================================

  # Risk-set-weighted covariate means (constant across rows)
  w_c <- ni / Lambda_C_est
  den <- sum(w_c)

  X_bar <- matrix(
    rep(colSums(XX * w_c) / den, each = nrow(X)),
    nrow = nrow(X)
  )

  if (sigma_sq_hat > 0) {
    B_bar <- matrix(
      rep(colSums(Bhat_subj * w_c) / den, each = nrow(X)),
      nrow = nrow(X)
    )
    design_M    <- cbind(X, Bhat)
    design_Mbar <- cbind(X_bar, B_bar)
  } else {
    design_M    <- X
    design_Mbar <- X_bar
  }

  # Filter to complete cases
  valid <- stats::complete.cases(design_M, y)

  H   <- design_M[valid, ] - design_Mbar[valid, ]
  M_v <- design_M[valid, ]
  y_v <- y[valid]

  LHS <- crossprod(H, M_v)
  RHS <- crossprod(H, y_v)
  beta_full <- tryCatch(
    as.vector(solve(LHS, RHS)),
    error = function(e) as.vector(MASS::ginv(LHS) %*% RHS)
  )

  p_X    <- ncol(X)
  p_full <- ncol(design_M)

  # ============================================
  # Step 6: Sandwich variance estimator
  # ============================================

  resid     <- as.vector(y_v - M_v %*% beta_full)
  scores    <- rowsum(H * resid, id_vec[valid])
  meat      <- crossprod(scores)
  bread_inv <- tryCatch(solve(LHS),
                        error = function(e) MASS::ginv(LHS))
  vcov_outcome <- bread_inv %*% meat %*% t(bread_inv)

  # ============================================
  # Unpack results
  # ============================================

  beta_hat <- beta_full[seq_len(p_X)]
  se_beta  <- sqrt(pmax(diag(vcov_outcome)[seq_len(p_X)], 0))
  names(beta_hat) <- long_covs
  names(se_beta)  <- long_covs

  if (sigma_sq_hat > 0 && p_full > p_X) {
    theta_hat <- beta_full[(p_X + 1L):p_full]
    se_theta  <- sqrt(pmax(diag(vcov_outcome)[(p_X + 1L):p_full], 0))
    theta_names <- colnames(Q0)
    if (!is.null(theta_names) && length(theta_names) == length(theta_hat)) {
      names(theta_hat) <- theta_names
      names(se_theta)  <- theta_names
    }
  } else {
    theta_hat <- NULL
    se_theta  <- NULL
  }

  vcov_names <- c(long_covs, if (sigma_sq_hat > 0 && p_full > p_X) colnames(Q0))
  if (length(vcov_names) == ncol(vcov_outcome))
    dimnames(vcov_outcome) <- list(vcov_names, vcov_names)

  # ============================================
  # Return
  # ============================================

  structure(list(
    gamma_hat    = gamma_hat,
    se_gamma     = se_gamma,
    vcov_gamma   = vcov_gamma,
    Lambda_C_est = Lambda_C_est,
    sigma_sq_hat = sigma_sq_hat,
    beta_hat     = beta_hat,
    se_beta      = se_beta,
    theta_hat    = theta_hat,
    se_theta     = se_theta,
    vcov_outcome = vcov_outcome,
    call         = call
  ), class = "Joint_modeling_visiting_and_longitudinal_Liang")
}

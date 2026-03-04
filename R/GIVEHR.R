#' Generalized Informative Visiting and Examination in EHR Data
#'
#' @description
#' Fits a three-stage semiparametric joint model for longitudinal electronic
#' health record (EHR) data that simultaneously accounts for an informative
#' visiting process and an informative observation process. Stage 1 estimates
#' the visiting-process intensity via partial likelihood with a log-normal
#' frailty. Stage 2 models the observation indicator with a probit regression
#' that may include latent-factor-linked random effects. Stage 3 estimates the
#' longitudinal outcome coefficients by weighted least squares with risk-set
#' centering. Sandwich variance estimators are provided for all three stages.
#'
#' @param data A \code{data.frame} in long (panel) format with one row per
#'   subject-time combination. Must contain columns for subject identifier,
#'   visit time, outcome, observation indicator, censoring time, and all
#'   specified covariates. Rows with \code{NA} in the time column are treated
#'   as unvisited time points.
#' @param id_col Character string giving the name of the subject identifier
#'   column in \code{data}. Default is \code{"id"}.
#' @param time_col Character string giving the name of the visit-time column
#'   in \code{data}. Default is \code{"time"}.
#' @param y_col Character string giving the name of the outcome column in
#'   \code{data}. Default is \code{"Y"}.
#' @param r_col Character string giving the name of the binary observation
#'   indicator column (\code{1} = outcome observed at visit, \code{0} =
#'   outcome missing at visit). Default is \code{"R"}.
#' @param censor_col Character string giving the name of the column containing
#'   each subject's censoring (end-of-follow-up) time. One value per subject
#'   is expected; the value on the first row for each subject is used. Default
#'   is \code{"C"}.
#' @param covars_visit_XV Character vector of covariate column names for the
#'   visiting-process intensity model. These enter the exponent of the
#'   multiplicative intensity. Default is \code{c("Z", "X")}.
#' @param covars_outcome_fixed_XY Character vector of covariate column names
#'   for the fixed-effect part of the longitudinal outcome model. No intercept
#'   is added automatically. Default is \code{c("Z", "X")}.
#' @param covars_outcome_link_ZY Character vector of covariate column names
#'   that interact with the conditional expectation of the latent factor
#'   \eqn{U_i} in the outcome model. An intercept column is prepended
#'   automatically. Including \code{time_col} in this vector adds a
#'   time-varying latent link. Default is \code{c("Z")}.
#' @param covars_obs_fixed_XO Character vector of covariate column names for
#'   the fixed-effect part of the observation-process probit model. An
#'   intercept is prepended automatically. Default is \code{c("Z", "X")}.
#' @param covars_obs_random_ZO Character vector of covariate column names that
#'   receive independent random effects in the observation process. An
#'   intercept is prepended automatically. Set to \code{NULL} or an empty
#'   vector to fit a fixed-effects-only probit. Default is \code{"Z"}.
#' @param time_precision Positive numeric scalar for time discretisation.
#'   Visit times are rounded to the nearest multiple of \code{time_precision}
#'   before constructing risk sets. Set to \code{NULL} or \code{0} to skip
#'   rounding. Default is \code{0.01}.
#'
#' @return A named list with three components, each itself a named list:
#'
#' \describe{
#'   \item{\code{visiting_process}}{A list with elements:
#'     \describe{
#'       \item{\code{gamma_V}}{Named numeric vector of visiting-process
#'         regression coefficients from partial likelihood.}
#'       \item{\code{se_gamma}}{Named numeric vector of sandwich standard
#'         errors for \code{gamma_V}.}
#'       \item{\code{vcov_gamma}}{Numeric matrix; the sandwich
#'         variance-covariance matrix for \code{gamma_V}.}
#'       \item{\code{Lambda0_hat}}{Named numeric vector of Breslow baseline
#'         hazard increments \eqn{d\hat\Lambda_0(t_k)} at each unique
#'         discretised event time.}
#'       \item{\code{sigma2_zeta}}{Numeric scalar; estimated log-normal
#'         frailty variance on the log scale.}
#'       \item{\code{mu_Ui}}{Named numeric vector of empirical Bayes
#'         posterior means for the standardised latent factor \eqn{U_i}.}
#'       \item{\code{sigma2_Ui}}{Named numeric vector of empirical Bayes
#'         posterior variances for \eqn{U_i}.}
#'       \item{\code{Ci}}{Named numeric vector of per-subject censoring
#'         times.}
#'       \item{\code{mi}}{Named integer vector of per-subject visit counts.}
#'     }
#'   }
#'   \item{\code{observation_process}}{A list with elements:
#'     \describe{
#'       \item{\code{alpha_fixed}}{Named numeric vector of fixed-effect
#'         probit coefficients (including intercept).}
#'       \item{\code{se_alpha}}{Named numeric vector of sandwich standard
#'         errors for \code{alpha_fixed}.}
#'       \item{\code{vcov_obs}}{Numeric matrix; sandwich
#'         variance-covariance matrix for all observation-process
#'         parameters.}
#'       \item{\code{delta_latent_link}}{(Present only when random effects
#'         are fitted.) Named numeric vector of latent-link coefficients
#'         \eqn{\delta}.}
#'       \item{\code{se_delta}}{(Present only when random effects are
#'         fitted.) Named numeric vector of sandwich standard errors for
#'         \code{delta_latent_link}.}
#'       \item{\code{Sigma_q}}{(Present only when random effects are
#'         fitted.) Named numeric vector of estimated random-effect
#'         standard deviations.}
#'     }
#'   }
#'   \item{\code{longitudinal_outcome}}{A list with elements:
#'     \describe{
#'       \item{\code{beta_Y}}{Named numeric vector of outcome fixed-effect
#'         coefficients.}
#'       \item{\code{theta_Y}}{Named numeric vector of outcome latent-link
#'         coefficients.}
#'       \item{\code{se_beta_naive}}{Named numeric vector of sandwich
#'         standard errors for \code{beta_Y}.}
#'       \item{\code{vcov_joint}}{Numeric matrix; sandwich
#'         variance-covariance matrix for the full coefficient vector
#'         \code{c(beta_Y, theta_Y)}.}
#'     }
#'   }
#' }
#'
#' @details
#' The GIVEHR estimator proceeds in three stages:
#'
#' \strong{Stage 1 - Visiting Process.}
#' The recurrent-event intensity for subject \eqn{i} is modelled as
#' \deqn{\lambda_i(t) = \eta_i \exp(\gamma^\top X_i^V)\,\lambda_0(t),}
#' where \eqn{\eta_i = \exp(\mu_0 + \sigma_\zeta U_i)} is a log-normal
#' frailty with \eqn{U_i \sim N(0,1)} and
#' \eqn{\mu_0 = -\sigma_\zeta^2 / 2} ensures \eqn{E[\eta_i] = 1}.
#' The coefficient vector \eqn{\gamma} is estimated by Newton-Raphson
#' on the partial likelihood score, the baseline hazard
#' \eqn{\lambda_0} by the Breslow estimator, and the frailty variance
#' \eqn{\sigma_\zeta^2} by a method-of-moments estimator. Subject-level
#' posterior summaries \eqn{(\hat\mu_{U_i}, \hat\sigma^2_{U_i})} are
#' obtained via Laplace approximation. A sandwich estimator
#' provides inference for \eqn{\gamma}.
#'
#' \strong{Stage 2 - Observation Process.}
#' At each visit the probability that the outcome is recorded is modelled as
#' \deqn{P(R_{ij}=1 \mid U_i, q_i) =
#'   \Phi\!\bigl(\alpha^\top X_{ij}^O +
#'   \delta^\top Z_{ij}^O \cdot U_i +
#'   q_i^\top Z_{ij}^O\bigr),}
#' where \eqn{q_i \sim N(0, \Sigma_q)} are independent random effects and
#' \eqn{\Phi} is the standard-normal CDF. The marginal log-likelihood
#' (integrated over \eqn{q_i} using the empirical Bayes posterior of
#' \eqn{U_i}) is maximised by \code{\link[stats]{nlminb}} with analytic
#' gradients. When the estimated frailty variance is negligible or no random
#' covariates are specified, a fixed-effects probit is fitted by iteratively
#' reweighted least squares. A clustered sandwich estimator provides
#' standard errors.
#'
#' \strong{Stage 3 - Longitudinal Outcome.}
#' The outcome model is
#' \deqn{Y_i(t) = \beta^\top X_i^Y(t) +
#'   \theta^\top Z_i^Y(t)\,\kappa_i(t) + \varepsilon_i(t),}
#' where \eqn{\kappa_i(t) = E[U_i \mid \text{observation history up to } t]}
#' is the conditional expectation of the latent factor given the
#' observation history. Estimation uses weighted least squares with
#' time-varying risk-set centering, where the weights combine the marginal
#' observation probability \eqn{\varpi_i(t)} and the visit intensity.
#' A sandwich estimator provides standard errors for \eqn{\beta}.
#'
#' @references
#' Yang, C.-H., Shi, X., and Mukherjee, B. (2026).
#' "Joint Modeling of Longitudinal EHR Data with Shared Random Effects for Informative Visiting and Observation Processes."
#' arXiv:2602.15374. doi:10.48550/arXiv.2602.15374.
#'
#' @seealso
#' \code{\link{GIVEHR_GH}} for a Gauss-Hermite quadrature variant with
#' online filtering;
#' \code{\link{GIVEHR_OU_PCL}} for the Ornstein-Uhlenbeck pairwise
#' composite likelihood extension;
#' \code{\link{EHRJoint}} for a simpler joint model without a separate
#' observation process;
#' \code{\link{sim_data_gen}} for simulating data compatible with this model.
#'
#' @examples
#' \dontrun{
#' # Simulate data
#' dat <- sim_data_gen(n = 500, seed = 42)
#'
#' # Fit the GIVEHR model
#' fit <- GIVEHR(
#'   data          = dat$long_data,
#'   id_col        = "id",
#'   time_col      = "time",
#'   y_col         = "Y",
#'   r_col         = "R",
#'   censor_col    = "C",
#'   covars_visit_XV         = c("Z", "X"),
#'   covars_outcome_fixed_XY = c("Z", "X"),
#'   covars_outcome_link_ZY  = "Z",
#'   covars_obs_fixed_XO     = c("Z", "X"),
#'   covars_obs_random_ZO    = "Z"
#' )
#'
#' # Inspect results
#' fit$visiting_process$gamma_V
#' fit$longitudinal_outcome$beta_Y
#' fit$longitudinal_outcome$se_beta_naive
#' summary(fit)
#' }
#'
#' @importFrom stats binomial coef dnorm glm nlminb plogis pnorm reformulate
#' @export
GIVEHR <- function(
    data,
    id_col                  = "id",
    time_col                = "time",
    y_col                   = "Y",
    r_col                   = "R",
    censor_col              = "C",
    covars_visit_XV         = c("Z", "X"),
    covars_outcome_fixed_XY = c("Z", "X"),
    covars_outcome_link_ZY  = c("Z"),
    covars_obs_fixed_XO     = c("Z", "X"),
    covars_obs_random_ZO    = "Z",
    time_precision          = 0.01
) {

  # ====================================================================
  # 0. Data prep
  # ====================================================================
  id_vec   <- data[[id_col]]
  time_vec <- data[[time_col]]
  y_vec    <- data[[y_col]]
  r_vec    <- data[[r_col]]
  n_row    <- length(id_vec)

  subj_first_idx <- !duplicated(id_vec)
  id_levels      <- id_vec[subj_first_idx]
  n_subj         <- length(id_levels)
  subj_index     <- match(id_vec, id_levels)
  sf             <- which(subj_first_idx)

  Ci_subj        <- data[[censor_col]][sf]
  visit_mask     <- !is.na(time_vec)
  visit_rows_idx <- which(visit_mask)
  n_visits_total <- length(visit_rows_idx)
  mi             <- tabulate(subj_index[visit_mask], nbins = n_subj)

  # Time discretization
  raw_times <- time_vec[visit_rows_idx]
  if (!is.null(time_precision) && time_precision > 0)
    raw_times <- round(raw_times / time_precision) * time_precision
  unique_times   <- sort(unique(raw_times))
  K              <- length(unique_times)
  event_time_idx <- match(raw_times, unique_times)
  event_subj_idx <- subj_index[visit_rows_idx]
  dN_counts      <- tabulate(event_time_idx, nbins = K)

  # Design matrices
  XV_subj <- as.matrix(data[sf, covars_visit_XV, drop = FALSE])
  p_V     <- ncol(XV_subj)

  build_dm <- function(vars, intercept = TRUE) {
    M <- as.matrix(data[, vars, drop = FALSE])
    if (intercept) M <- cbind(`(Intercept)` = 1, M)
    M
  }
  XO_all <- build_dm(covars_obs_fixed_XO, intercept = TRUE); p_XO <- ncol(XO_all)
  has_random <- !is.null(covars_obs_random_ZO) && length(covars_obs_random_ZO) > 0
  if (has_random) { ZO_all <- build_dm(covars_obs_random_ZO, intercept = TRUE); q_dim <- ncol(ZO_all)
  } else { ZO_all <- NULL; q_dim <- 0L }
  XY_all <- build_dm(covars_outcome_fixed_XY, intercept = FALSE); p_XY <- ncol(XY_all)

  # Helper: tcrossprod(X, beta_mat) replaces X %*% beta — 1.5x faster at scale
  # beta_mat must be matrix(beta, nrow=1)
  mat_vec <- function(X, beta) drop(tcrossprod(X, matrix(beta, nrow = 1)))

  # ====================================================================
  # Risk-set engine
  # ====================================================================
  if (!is.null(time_precision) && time_precision > 0) {
    Ci_disc <- floor(Ci_subj / time_precision) * time_precision
  } else Ci_disc <- Ci_subj
  last_risk_idx <- findInterval(Ci_disc, unique_times)
  valid_subj    <- which(last_risk_idx >= 1L)
  lri_valid     <- last_risk_idx[valid_subj]
  n_valid       <- length(valid_subj)

  rs_accum_scalar <- function(w) {
    acc <- numeric(K)
    if (n_valid > 0L) {
      rs <- rowsum(w[valid_subj], lri_valid, reorder = FALSE)
      acc[as.integer(rownames(rs))] <- rs[, 1L]
    }
    rev(cumsum(rev(acc)))
  }

  rs_accum_matrix <- function(w, M) {
    p <- ncol(M); acc <- matrix(0, K, p)
    if (n_valid > 0L) {
      rs <- rowsum(w[valid_subj] * M[valid_subj, , drop = FALSE],
                   lri_valid, reorder = FALSE)
      acc[as.integer(rownames(rs)), ] <- as.matrix(rs)
    }
    for (j in seq_len(p)) acc[, j] <- rev(cumsum(rev(acc[, j])))
    acc
  }

  rs_accum_01 <- function(w, M) {
    p <- ncol(M)
    aug <- cbind(w[valid_subj], w[valid_subj] * M[valid_subj, , drop = FALSE])
    acc <- matrix(0, K, 1L + p)
    if (n_valid > 0L) {
      rs <- rowsum(aug, lri_valid, reorder = FALSE)
      acc[as.integer(rownames(rs)), ] <- as.matrix(rs)
    }
    for (j in seq_len(1L + p)) acc[, j] <- rev(cumsum(rev(acc[, j])))
    list(s0 = acc[, 1L], sM = acc[, -1L, drop = FALSE])
  }

  # ====================================================================
  # Stage 1: AG Partial Likelihood
  # ====================================================================
  score_hess <- function(gamma) {
    rs_vec <- exp(mat_vec(XV_subj, gamma))
    sh <- rs_accum_01(rs_vec, XV_subj)
    S0 <- sh$s0; S1 <- sh$sM
    inv_S0 <- 1 / pmax(S0, 1e-30)
    Xbar <- S1 * inv_S0

    Xe     <- XV_subj[event_subj_idx, , drop = FALSE]
    Xbar_e <- Xbar[event_time_idx, , drop = FALSE]
    sc     <- drop(crossprod(rep(1, n_visits_total), Xe - Xbar_e))

    H <- matrix(0, p_V, p_V)
    for (a in seq_len(p_V)) {
      S2_a <- rs_accum_matrix(rs_vec * XV_subj[, a], XV_subj)
      for (b in a:p_V) {
        val <- -sum(dN_counts * (S2_a[, b] * inv_S0 - Xbar[, a] * Xbar[, b]))
        H[a, b] <- val; H[b, a] <- val
      }
    }
    list(score = sc, hessian = H, S0 = S0, Xbar = Xbar, rs_vec = rs_vec)
  }

  gamma_hat <- numeric(p_V)
  sh_final  <- NULL
  for (it in 1:50) {
    sh_final <- score_hess(gamma_hat)
    if (max(abs(sh_final$score)) < 1e-8) break
    step <- tryCatch(solve(sh_final$hessian, sh_final$score), error = function(e) numeric(p_V))
    step[!is.finite(step)] <- 0
    step <- pmax(pmin(step, 3), -3)
    gamma_hat <- gamma_hat - step
  }
  names(gamma_hat) <- colnames(XV_subj)

  # -- Stage 1 Sandwich Variance (Lin-Wei-Ying) --
  # Bread = -Hessian at final iteration
  Info_gamma <- -sh_final$hessian
  inv_Info   <- tryCatch(solve(Info_gamma), error = function(e) base::qr.solve(Info_gamma))

  # Breslow
  rs_final <- sh_final$rs_vec
  S0_final <- sh_final$S0
  Xbar_final <- sh_final$Xbar
  dLam0    <- dN_counts / pmax(S0_final, 1e-30)
  cumLam0  <- cumsum(dLam0)

  Lam0_Ci <- numeric(n_subj)
  Lam0_Ci[valid_subj] <- cumLam0[lri_valid]

  # Per-subject score residual L_i (p_V-vector each):
  #   L_i = [event term] - rs_i * [compensator term]
  # event_i = sum_{j: visits of i} (X_i - Xbar(t_ij))
  # comp_i  = sum_{k: t_k <= C_i} (X_i - Xbar(t_k)) * dLambda0(t_k)
  #         = X_i * Lambda0(C_i) - A(C_i)
  # where A(t) = cumsum_k [Xbar(t_k) * dLambda0(t_k)]

  # Event term: aggregate (X_i - Xbar_event) by subject
  resid_event <- XV_subj[event_subj_idx, , drop = FALSE] -
    Xbar_final[event_time_idx, , drop = FALSE]                    # n_visits x p_V
  event_term_i <- matrix(0, n_subj, p_V)
  if (n_visits_total > 0L) {
    event_term_i <- rowsum(resid_event, event_subj_idx, reorder = FALSE)
    # rowsum returns only subjects that have events; pad to n_subj
    evt_ids  <- as.integer(rownames(event_term_i))
    tmp      <- matrix(0, n_subj, p_V)
    tmp[evt_ids, ] <- event_term_i
    event_term_i   <- tmp
  }

  # Compensator term: A(C_i) = cumsum(Xbar * dLambda0) evaluated at last_risk_idx
  Xbar_dLam <- Xbar_final * dLam0                                 # K x p_V
  cumA <- apply(Xbar_dLam, 2, cumsum)                              # K x p_V
  A_Ci <- matrix(0, n_subj, p_V)
  A_Ci[valid_subj, ] <- cumA[lri_valid, , drop = FALSE]

  comp_term_i <- XV_subj * Lam0_Ci - A_Ci                         # n_subj x p_V

  L_i <- event_term_i - rs_final * comp_term_i                    # n_subj x p_V
  Meat_gamma <- crossprod(L_i)                                     # p_V x p_V
  vcov_gamma <- inv_Info %*% Meat_gamma %*% inv_Info
  se_gamma   <- sqrt(pmax(diag(vcov_gamma), 0))

  # -- Continue Stage 1: frailty + EB --
  nu_i  <- rs_final * Lam0_Ci
  mm_n  <- sum(mi * (mi - 1L)) - sum(nu_i^2)
  mm_d  <- sum(nu_i^2)
  var_eta <- if (mm_d > 0) max(mm_n / mm_d, 0) else 0
  s2z     <- if (var_eta > 0) log1p(var_eta) else 0
  mu0z    <- -0.5 * s2z

  if (s2z > 0) {
    nu_safe <- pmax(nu_i, 1e-12); log_nu <- log(nu_safe)
    mode_v <- ifelse(mi == 0L, mu0z - 0.5, log((mi + 0.5) / nu_safe))
    for (k in 1:60) {
      e_term <- exp(log_nu + mode_v)
      g <- mi - e_term - (mode_v - mu0z) / s2z
      Hv <- -e_term - 1 / s2z
      step <- g / Hv; step[!is.finite(step)] <- 0
      step <- pmax(pmin(step, 2), -2)
      mode_new <- mode_v - step
      if (max(abs(step)) < 1e-8) break
      mode_v <- mode_new
    }
    ne_mode <- exp(log_nu + mode_v)
    var_le  <- 1 / (ne_mode + 1 / s2z)
    sq_s2z  <- sqrt(s2z)
    mu_U    <- (mode_v - mu0z) / sq_s2z
    s2_U    <- var_le / s2z
  } else {
    mu_U <- rep(0, n_subj); s2_U <- rep(0, n_subj)
  }
  mu_U_r <- mu_U[subj_index]; s2_U_r <- s2_U[subj_index]

  # ====================================================================
  # Stage 2: Observation Process
  # ====================================================================
  obs_idx <- which(!is.na(r_vec) & visit_mask)
  n_obs   <- length(obs_idx)
  XO_obs  <- XO_all[obs_idx, , drop = FALSE]
  R_obs   <- r_vec[obs_idx]
  is_fz   <- (var_eta < 1e-6)
  do_rand <- has_random && !is_fz

  if (do_rand) {
    ZO_obs    <- ZO_all[obs_idx, , drop = FALSE]
    ZO_obs_sq <- ZO_obs^2
    mu_U_obs  <- mu_U_r[obs_idx]
    s2_U_obs  <- s2_U_r[obs_idx]
    sign_R    <- 2 * R_obs - 1
  }

  if (!do_rand) {
    alpha_hat <- numeric(p_XO)
    for (.iwls in 1:15) {
      eta <- mat_vec(XO_obs, alpha_hat)
      eta[eta > 5] <- 5; eta[eta < -5] <- -5
      mu <- pnorm(eta); mu[mu < 1e-6] <- 1e-6; mu[mu > 1-1e-6] <- 1-1e-6
      dphi <- dnorm(eta); w <- dphi^2 / (mu * (1 - mu))
      z <- eta + (R_obs - mu) / dphi; WX <- w * XO_obs
      a_new <- drop(solve(crossprod(XO_obs, WX), crossprod(WX, z)))
      if (max(abs(a_new - alpha_hat)) < 1e-8) { alpha_hat <- a_new; break }
      alpha_hat <- a_new
    }
    names(alpha_hat) <- colnames(XO_all)
    delta_hat <- setNames(numeric(0), character(0))
    Sq_hat    <- setNames(numeric(0), character(0))

    # Stage 2 vcov (probit, no random effects): standard sandwich
    eta_f <- mat_vec(XO_obs, alpha_hat)
    eta_f[eta_f > 5] <- 5; eta_f[eta_f < -5] <- -5
    mu_f  <- pnorm(eta_f); mu_f[mu_f < 1e-6] <- 1e-6; mu_f[mu_f > 1-1e-6] <- 1-1e-6
    dphi_f <- dnorm(eta_f)
    w_f    <- dphi_f^2 / (mu_f * (1 - mu_f))
    # Info = X' W X
    Info_obs <- crossprod(XO_obs, w_f * XO_obs)
    inv_Info_obs <- tryCatch(solve(Info_obs), error = function(e) base::qr.solve(Info_obs))
    # Per-obs score: (R - mu) * dphi/(mu*(1-mu)) * X
    s_obs_row <- ((R_obs - mu_f) * dphi_f / (mu_f * (1 - mu_f))) * XO_obs
    s_obs_i   <- rowsum(s_obs_row, subj_index[obs_idx])
    Meat_obs  <- crossprod(s_obs_i)
    vcov_obs_all <- inv_Info_obs %*% Meat_obs %*% inv_Info_obs

  } else {
    # Full GLM init (matches original — glm.fit is not the bottleneck)
    fit0 <- stats::glm(stats::reformulate(covars_obs_fixed_XO, response = r_col),
                       data = data[obs_idx, , drop = FALSE],
                       family = stats::binomial(link = "probit"))
    alpha_init <- numeric(p_XO); names(alpha_init) <- colnames(XO_obs)
    tmp <- stats::coef(fit0); alpha_init[names(tmp)] <- tmp
    n_par <- p_XO + 2L * q_dim

    # fn_and_gr: returns nll with gradient attribute
    # Also optionally returns per-obs score matrix (for vcov)
    fn_and_gr <- function(par, return_scores = FALSE) {
      alpha  <- par[1:p_XO]
      delta  <- par[(p_XO + 1):(p_XO + q_dim)]
      log_sq <- par[(p_XO + q_dim + 1):n_par]
      sq     <- exp(log_sq)

      aX  <- mat_vec(XO_obs, alpha)
      dZ  <- mat_vec(ZO_obs, delta)
      vrs <- mat_vec(ZO_obs_sq, sq)

      tv  <- 1 + vrs + dZ * dZ * s2_U_obs
      isd <- 1 / sqrt(tv)
      k   <- (aX + dZ * mu_U_obs) * isd
      k[k > 8] <- 8; k[k < -8] <- -8

      sk      <- sign_R * k
      log_Phi <- pnorm(sk, log.p = TRUE)
      nll     <- -sum(log_Phi)

      mills <- exp(dnorm(sk, log = TRUE) - log_Phi) * sign_R
      k_over_tv <- k / tv

      # Per-obs gradient components (n_obs vectors)
      mills_isd <- mills * isd
      coef2     <- mills * k_over_tv * s2_U_obs * dZ
      coef1     <- mills * mu_U_obs * isd
      coef_lsq  <- mills * (-0.5 * k_over_tv)

      if (!return_scores) {
        g_alpha <- -drop(crossprod(mills_isd, XO_obs))
        g_delta <- -drop(crossprod(coef1, ZO_obs)) + drop(crossprod(coef2, ZO_obs))
        g_logsq <- -sq * drop(crossprod(coef_lsq, ZO_obs_sq))
        attr(nll, "gradient") <- c(g_alpha, g_delta, g_logsq)
        return(nll)
      }

      # Return per-observation score matrix (n_obs x n_par) for sandwich
      s_alpha <- -mills_isd * XO_obs                                    # n_obs x p_XO
      s_delta <- -(coef1 * ZO_obs) + (coef2 * ZO_obs)                  # n_obs x q_dim
      s_logsq <- -(coef_lsq * ZO_obs_sq) * rep(sq, each = n_obs)       # n_obs x q_dim
      score_mat <- cbind(s_alpha, s_delta, s_logsq)

      attr(nll, "gradient") <- drop(crossprod(rep(1, n_obs), score_mat))
      attr(nll, "scores")   <- score_mat
      nll
    }

    .cp <- NULL; .cv <- NULL
    obj_fn <- function(par) {
      if (!identical(par, .cp)) { .cv <<- fn_and_gr(par); .cp <<- par }
      as.numeric(.cv)
    }
    gr_fn <- function(par) {
      if (!identical(par, .cp)) { .cv <<- fn_and_gr(par); .cp <<- par }
      attr(.cv, "gradient")
    }

    par0 <- c(alpha_init, rep(0, q_dim), rep(log(0.2), q_dim))
    opt  <- stats::nlminb(par0, obj_fn, gr_fn, control = list(iter.max = 200))

    alpha_hat <- setNames(opt$par[1:p_XO], colnames(XO_obs))
    delta_hat <- setNames(opt$par[(p_XO + 1):(p_XO + q_dim)], colnames(ZO_obs))
    Sq_hat    <- setNames(exp(opt$par[(p_XO + q_dim + 1):n_par]), colnames(ZO_obs))

    # -- Stage 2 Sandwich Variance --
    # Meat: per-subject clustered scores
    final_eval <- fn_and_gr(opt$par, return_scores = TRUE)
    score_obs_mat <- attr(final_eval, "scores")                     # n_obs x n_par
    score_subj_i  <- rowsum(score_obs_mat, subj_index[obs_idx])     # n_unique_subj x n_par
    Meat_obs      <- crossprod(score_subj_i)                        # n_par x n_par

    # Bread: forward-difference Hessian (n_par calls, not 2*n_par)
    gr_at_opt <- attr(final_eval, "gradient")
    eps_h     <- 1e-5
    H_obs     <- matrix(0, n_par, n_par)
    for (j in seq_len(n_par)) {
      e_j <- numeric(n_par); e_j[j] <- eps_h
      gr_plus    <- attr(fn_and_gr(opt$par + e_j), "gradient")
      H_obs[, j] <- (gr_plus - gr_at_opt) / eps_h
    }
    H_obs <- 0.5 * (H_obs + t(H_obs))
    inv_H_obs <- tryCatch(solve(H_obs), error = function(e) base::qr.solve(H_obs))
    vcov_obs_all <- inv_H_obs %*% Meat_obs %*% inv_H_obs
  }

  # ====================================================================
  # Compute varpi and kappa (all rows)
  # ====================================================================
  aX_all <- mat_vec(XO_all, alpha_hat)
  if (do_rand) {
    dZ_all  <- mat_vec(ZO_all, delta_hat)
    vrs_all <- mat_vec(ZO_all^2, Sq_hat)
  } else {
    dZ_all <- numeric(n_row); vrs_all <- numeric(n_row)
  }
  tv_all <- 1 + vrs_all + dZ_all * dZ_all * s2_U_r
  tv_all[tv_all < 1 + 1e-12] <- 1 + 1e-12
  isd_all <- 1 / sqrt(tv_all)
  k_all   <- (aX_all + dZ_all * mu_U_r) * isd_all
  k_all[k_all > 12] <- 12; k_all[k_all < -12] <- -12

  varpi <- pnorm(k_all); phi_k <- dnorm(k_all)
  varpi_safe <- varpi; varpi_safe[varpi_safe < 1e-12] <- 1e-12
  kappa <- mu_U_r + dZ_all * s2_U_r * isd_all * phi_k / varpi_safe

  # ====================================================================
  # Stage 3: WLS with time-varying centering
  # ====================================================================
  use_time <- time_col %in% covars_outcome_link_ZY
  zy_vars  <- setdiff(covars_outcome_link_ZY, time_col)
  if (length(zy_vars) > 0) {
    ZY_row <- cbind(1, as.matrix(data[, zy_vars, drop = FALSE]))
    colnames(ZY_row)[1] <- "Intercept"
  } else {
    ZY_row <- matrix(1, n_row, 1, dimnames = list(NULL, "Intercept"))
  }
  if (use_time) ZY_row <- cbind(ZY_row, time = time_vec)

  B_all <- kappa * ZY_row

  varpi_s <- varpi[sf]; kappa_s <- kappa[sf]
  w_s <- varpi_s * mi / pmax(Lam0_Ci, 1e-30)
  w_s[mi == 0L | Lam0_Ci <= 0] <- 0

  XY_s <- XY_all[sf, , drop = FALSE]
  ZY_s <- ZY_row[sf, , drop = FALSE]
  ZY_s_notime <- if (use_time) ZY_s[, -ncol(ZY_s), drop = FALSE] else ZY_s
  B_s  <- kappa_s * ZY_s_notime; p_B <- ncol(B_s)

  dw_k    <- rs_accum_scalar(w_s)
  inv_dw  <- 1 / pmax(dw_k, 1e-30)
  XYbar_k <- rs_accum_matrix(w_s, XY_s) * inv_dw
  Bbar_k  <- rs_accum_matrix(w_s, B_s)  * inv_dw
  if (use_time) kbar_k <- rs_accum_scalar(w_s * kappa_s) * inv_dw

  row_tidx <- match(raw_times, unique_times)

  eo_mask <- visit_mask & (r_vec == 1L) & !is.na(y_vec)
  eo_idx  <- which(eo_mask); n_eo <- length(eo_idx)
  if (n_eo == 0L) stop("No observed outcome rows.")

  eo_in_visits <- match(eo_idx, visit_rows_idx)
  eo_time_k    <- event_time_idx[eo_in_visits]

  D_eo <- cbind(XY_all[eo_idx, , drop = FALSE], B_all[eo_idx, , drop = FALSE])
  XYbar_eo <- XYbar_k[eo_time_k, , drop = FALSE]
  Bbar_eo  <- Bbar_k[eo_time_k, , drop = FALSE]
  if (use_time) {
    Bbar_eo <- cbind(Bbar_eo, kbar_k[eo_time_k] * time_vec[eo_idx])
  }
  Dbar_eo <- cbind(XYbar_eo, Bbar_eo)

  H_eo <- D_eo - Dbar_eo; Y_eo <- y_vec[eo_idx]
  LHS <- crossprod(H_eo, D_eo); RHS <- crossprod(H_eo, Y_eo)
  coefs <- drop(tryCatch(solve(LHS, RHS), error = function(e) base::qr.solve(LHS, RHS)))

  beta_hat  <- coefs[1:p_XY]; theta_hat <- coefs[-(1:p_XY)]
  names(beta_hat)  <- colnames(XY_all)
  names(theta_hat) <- colnames(B_all)

  # Stage 3 sandwich
  res_eo  <- Y_eo - drop(tcrossprod(D_eo, matrix(coefs, nrow = 1)))
  score_i <- rowsum(H_eo * res_eo, subj_index[eo_idx])
  Meat3   <- crossprod(score_i)
  invB3   <- tryCatch(solve(LHS), error = function(e) base::qr.solve(LHS))
  V3      <- invB3 %*% Meat3 %*% t(invB3)
  se_beta <- sqrt(pmax(diag(V3)[1:p_XY], 0))

  # ====================================================================
  # Return
  # ====================================================================
  # Stage 2 vcov: extract SE for each parameter group
  par_names_obs <- c(names(alpha_hat),
                     if (do_rand) paste0("delta.", names(delta_hat)),
                     if (do_rand) paste0("logSq.", names(Sq_hat)))
  if (exists("vcov_obs_all") && nrow(vcov_obs_all) > 0)
    dimnames(vcov_obs_all) <- list(par_names_obs, par_names_obs)

  se_alpha <- sqrt(pmax(diag(vcov_obs_all)[1:p_XO], 0))
  names(se_alpha) <- names(alpha_hat)

  obs_out <- list(
    alpha_fixed = alpha_hat,
    se_alpha    = se_alpha,
    vcov_obs    = vcov_obs_all
  )
  if (do_rand) {
    obs_out$delta_latent_link <- delta_hat
    obs_out$Sigma_q           <- Sq_hat
    se_delta <- sqrt(pmax(diag(vcov_obs_all)[(p_XO + 1):(p_XO + q_dim)], 0))
    names(se_delta) <- names(delta_hat)
    obs_out$se_delta <- se_delta
  }

  list(
    visiting_process = list(
      gamma_V     = gamma_hat,
      se_gamma    = setNames(se_gamma, names(gamma_hat)),
      vcov_gamma  = vcov_gamma,
      Lambda0_hat = setNames(dLam0, as.character(unique_times)),
      sigma2_zeta = s2z,
      mu_Ui       = setNames(mu_U, id_levels),
      sigma2_Ui   = setNames(s2_U, id_levels),
      Ci          = setNames(Ci_subj, id_levels),
      mi          = setNames(mi, id_levels)
    ),
    observation_process = obs_out,
    longitudinal_outcome = list(
      beta_Y        = beta_hat,
      theta_Y       = theta_hat,
      se_beta_naive = setNames(se_beta, names(beta_hat)),
      vcov_joint    = V3
    )
  )
}

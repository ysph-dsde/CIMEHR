#' Clinical Informative Missingness for Electronic Health Records (CIMEHR)
#' with Ornstein--Uhlenbeck (OU) Pairwise Composite Likelihood (PCL) Extension
#' -- Time-Varying Variant
#'
#' @description
#' Fits a three-stage semiparametric joint model for longitudinal Electronic
#' Health Record (EHR) data, extending \code{\link{CIMEHR}} by modelling
#' temporal serial correlation in the observation process with a stationary
#' Ornstein--Uhlenbeck (OU) process and estimating the observation-process
#' parameters via Pairwise Composite Likelihood (PCL).
#'
#' Stage 1 estimates the visiting-process intensity via Andersen--Gill
#' partial likelihood with a log-normal frailty (identical to
#' \code{\link{CIMEHR}}). Stage 2 maximises a pairwise composite
#' log-likelihood over adjacent (or all) within-subject visit pairs,
#' integrating over the bivariate probit structure induced by the OU
#' process. An analytic gradient is provided for efficiency. Stage 3
#' estimates the longitudinal outcome coefficients by weighted least
#' squares with risk-set centering and a robust sandwich variance
#' estimator.
#'
#' @param data A \code{data.frame} in long (panel) format with one row per
#'   subject-time combination. Must contain columns for subject identifier,
#'   visit time, outcome, observation indicator, and all specified
#'   covariates. Rows with \code{NA} in the time column are treated as
#'   unvisited time points.
#' @param id_col Character string giving the name of the subject identifier
#'   column. Default is \code{"id"}.
#' @param time_col Character string giving the name of the visit-time
#'   column. Default is \code{"time"}.
#' @param y_col Character string giving the name of the outcome column.
#'   Default is \code{"Y"}.
#' @param r_col Character string giving the name of the binary observation
#'   indicator column (\code{1} = outcome observed, \code{0} = outcome
#'   missing at visit). Default is \code{"R"}.
#' @param covars_visit_XV Character vector of covariate column names for the
#'   visiting-process intensity model. Default is \code{c("Z", "X")}.
#' @param covars_outcome_fixed_XY Character vector of covariate column names
#'   for the fixed-effect part of the longitudinal outcome model. No
#'   intercept is added automatically. Default is \code{c("Z", "X")}.
#' @param covars_outcome_random_link_ZY Character vector of covariate column names
#'   that interact with the conditional expectation of the latent factor
#'   \eqn{U_i} in the outcome model. An intercept is prepended
#'   automatically. Including \code{time_col} adds a time-varying link.
#'   Default is \code{c("Z")}.
#' @param covars_obs_fixed_XO Character vector of covariate column names for
#'   the fixed-effect part of the observation-process probit model. An
#'   intercept is prepended automatically. Default is \code{c("Z", "X")}.
#' @param covars_obs_random_link_ZO Character vector of covariate column names
#'   receiving independent random effects in the observation process. An
#'   intercept is prepended automatically. Set to \code{NULL} for a
#'   fixed-effects-only probit (in which case PCL is skipped and a
#'   standard probit generalized linear model (GLM) is fitted). Default is \code{"Z"}.
#' @param pair_type Character string specifying which within-subject visit
#'   pairs to use for the pairwise composite likelihood. One of:
#'   \describe{
#'     \item{\code{"adjacent"}}{(Default.) Use only temporally adjacent
#'       pairs \eqn{(j, j+1)}. This is computationally cheaper and
#'       typically sufficient.}
#'     \item{\code{"all"}}{Use all \eqn{\binom{m_i}{2}} pairs per
#'       subject. More informative but \eqn{O(m_i^2)} per subject.}
#'   }
#' @param phi_init Numeric scalar giving the starting value for the OU
#'   decay parameter \eqn{\phi} in the optimisation. Default is
#'   \code{0.1}.
#' @param phi_bounds Numeric vector of length 2 giving the lower and upper
#'   bounds for \eqn{\phi} (on the original scale; internally
#'   log-transformed). Default is \code{c(1e-4, 100)}.
#' @param optim_method Character string specifying the optimisation method
#'   passed to \code{\link[stats]{optim}}. Default is \code{"L-BFGS-B"}.
#' @param verbose Logical; if \code{TRUE}, progress messages are printed
#'   during estimation. Default is \code{FALSE}.
#'
#' @return An object of class \code{"CIMEHR_timevarying_ou"}, which is a named list
#' with the following components:
#'
#' \describe{
#'   \item{\code{visiting_process}}{A list with elements:
#'     \describe{
#'       \item{\code{gamma_V}}{Named numeric vector of visiting-process
#'         regression coefficients.}
#'       \item{\code{Lambda0_hat}}{Named numeric vector of Breslow baseline
#'         hazard increments at each unique event time.}
#'       \item{\code{sigma2_zeta}}{Numeric scalar; estimated log-normal
#'         frailty variance on the log scale.}
#'       \item{\code{mu_Ui}}{Named numeric vector of empirical Bayes
#'         posterior means for \eqn{U_i}.}
#'       \item{\code{sigma2_Ui}}{Named numeric vector of empirical Bayes
#'         posterior variances for \eqn{U_i}.}
#'     }
#'   }
#'   \item{\code{observation_process}}{A list with elements:
#'     \describe{
#'       \item{\code{alpha_fixed}}{Named numeric vector of fixed-effect
#'         probit coefficients (including intercept).}
#'       \item{\code{delta_latent_random_link}}{(Present only when random effects
#'         are fitted.) Named numeric vector of latent-link coefficients
#'         \eqn{\delta}.}
#'       \item{\code{Sigma_q}}{(Present only when random effects are
#'         fitted.) Named numeric vector of estimated random-effect
#'         variances \eqn{\sigma_{q_l}^2}.}
#'       \item{\code{phi_OU}}{(Present only when random effects are
#'         fitted.) Numeric scalar; estimated OU decay parameter
#'         \eqn{\phi}. Larger values indicate faster decay of serial
#'         correlation.}
#'       \item{\code{pcl_value}}{(Present only when random effects are
#'         fitted.) Numeric scalar; minimised negative pairwise composite
#'         log-likelihood.}
#'       \item{\code{convergence}}{(Present only when random effects are
#'         fitted.) Integer convergence code from
#'         \code{\link[stats]{optim}}.}
#'     }
#'   }
#'   \item{\code{longitudinal_outcome}}{A list with elements:
#'     \describe{
#'       \item{\code{beta_Y}}{Named numeric vector of outcome fixed-effect
#'         coefficients.}
#'       \item{\code{theta_Y}}{Named numeric vector of outcome latent-link
#'         coefficients.}
#'       \item{\code{se_beta}}{Named numeric vector of sandwich standard
#'         errors for \code{beta_Y}.}
#'       \item{\code{vcov_robust}}{Numeric matrix; sandwich
#'         variance--covariance matrix for
#'         \code{c(beta_Y, theta_Y)}.}
#'     }
#'   }
#'   \item{\code{model_info}}{A list with elements:
#'     \describe{
#'       \item{\code{n_subjects}}{Integer; number of unique subjects.}
#'       \item{\code{n_visits}}{Integer; total number of visits.}
#'       \item{\code{n_observed}}{Integer; total number of visits with
#'         \code{R = 1}.}
#'       \item{\code{pair_type}}{Character string echoing the
#'         \code{pair_type} argument.}
#'       \item{\code{median_gap}}{Numeric scalar; median time gap between
#'         consecutive visits within subjects.}
#'     }
#'   }
#'   \item{\code{call}}{The matched call.}
#' }
#'
#' @details
#' The model extends \code{\link{CIMEHR}} by introducing an
#' Ornstein--Uhlenbeck (OU) process \eqn{\xi_i(t)} in the observation model:
#' \deqn{\tilde R_{ij} = \alpha^\top X_{ij}^O + \delta^\top Z_{ij}^O
#'   \cdot U_i + q_i^\top Z_{ij}^O + \xi_i(t_{ij}),}
#' where \eqn{R_{ij} = I(\tilde R_{ij} > 0)}, \eqn{q_i \sim N(0,
#' \Sigma_q)}, and \eqn{\xi_i} is a stationary OU process with unit
#' marginal variance and autocorrelation
#' \eqn{\mathrm{Cor}(\xi_i(s), \xi_i(t)) = \exp(-\phi |t - s|)}. The
#' OU component captures temporal serial correlation in the observation
#' process that is not explained by covariates or the shared latent
#' factor.
#'
#' \strong{Stage 2 -- Pairwise Composite Likelihood.}
#' Because the full likelihood involves high-dimensional integrals over
#' \eqn{\xi_i}, estimation proceeds via pairwise composite likelihood.
#' For each within-subject visit pair \eqn{(j, k)} the bivariate probit
#' probability \eqn{P(R_{ij} = r_j, R_{ik} = r_k)} is computed using
#' \code{\link[pbivnorm]{pbivnorm}}, and the composite log-likelihood is
#' \deqn{cl(\theta) = \sum_i \sum_{j < k} \log \Phi_2(d_j z_j,\, d_k
#'   z_k;\; d_j d_k \rho_{jk}),}
#' where \eqn{d_j = 2 R_{ij} - 1}, \eqn{z_j} is the standardised
#' marginal probit argument, and \eqn{\rho_{jk}} combines the latent
#' factor, random-effect, and OU contributions to the pairwise
#' correlation. An analytic gradient is provided and the composite
#' log-likelihood is maximised with \code{\link[stats]{optim}} using the
#' specified \code{optim_method}.
#'
#' \strong{Stages 1 and 3} are identical to \code{\link{CIMEHR}}, except
#' that Stage 3 uses a static (non-filtered) version of \eqn{\omega_i(t)}
#' and \eqn{\kappa_i(t)}.
#'
#' S3 methods \code{print}, \code{summary}, \code{coef}, \code{vcov}, and
#' \code{confint} are available for objects of class
#' \code{"CIMEHR_timevarying_ou"}.
#'
#' @references
#' Yang, C.-H., Shi, X., and Mukherjee, B. (2026).
#' "Joint Modeling of Longitudinal EHR Data with Shared Random Effects for Informative Visiting and Observation Processes."
#' arXiv:2602.15374. \doi{10.48550/arXiv.2602.15374}.
#'
#'
#' @seealso
#' \code{\link{CIMEHR}} for the base version without Ornstein--Uhlenbeck (OU) serial
#' correlation;
#' \code{\link{CIMEHR_timevarying_ghq}} for the Gauss--Hermite quadrature variant;
#' \code{\link{sim_data_gen}} for simulating data compatible with this
#' model;
#'
#' @examples
#' dat <- sim_ehr_data[sim_ehr_data$id <= 80, ]
#' fit <- CIMEHR_timevarying_ou(data = dat, pair_type = "adjacent")
#' print(fit)
#'
#' @importFrom stats as.formula binomial coef dnorm glm model.matrix
#'   optim pnorm qnorm reformulate setNames median
#' @importFrom utils combn
#' @export
CIMEHR_timevarying_ou <- function(
    data,
    id_col                  = "id",
    time_col                = "time",
    y_col                   = "Y",
    r_col                   = "R",
    covars_visit_XV         = c("Z", "X"),
    covars_outcome_fixed_XY = c("Z", "X"),
    covars_outcome_random_link_ZY  = c("Z"),
    covars_obs_fixed_XO     = c("Z", "X"),
    covars_obs_random_link_ZO    = "Z",
    pair_type               = c("adjacent", "all"),
    phi_init                = 0.1,
    phi_bounds              = c(1e-4, 100),
    optim_method            = "L-BFGS-B",
    verbose                 = FALSE
) {

  # Capture call
  call <- match.call()
  pair_type <- match.arg(pair_type)

  # ===========================================================================
  # Input Validation
  # ===========================================================================
  if (!is.data.frame(data)) {
    stop("`data` must be a data.frame.", call. = FALSE)
  }

  required_cols <- c(id_col, time_col, y_col, r_col)
  missing_cols <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0) {
    stop(sprintf("Required columns not found: %s",
                 paste(missing_cols, collapse = ", ")), call. = FALSE)
  }

  all_covars <- unique(c(
    covars_visit_XV, covars_outcome_fixed_XY, covars_outcome_random_link_ZY,
    covars_obs_fixed_XO, covars_obs_random_link_ZO
  ))
  all_covars <- setdiff(all_covars, time_col)
  missing_covars <- setdiff(all_covars, names(data))
  if (length(missing_covars) > 0) {
    stop(sprintf("Covariate columns not found: %s",
                 paste(missing_covars, collapse = ", ")), call. = FALSE)
  }

  # ===========================================================================
  # Load pbivnorm
  # ===========================================================================
  if (!requireNamespace("pbivnorm", quietly = TRUE)) {
    stop("Package 'pbivnorm' is required. Install it with: install.packages('pbivnorm')",
         call. = FALSE)
  }
  pbvnorm_fast <- pbivnorm::pbivnorm

  # ===========================================================================
  # 0. Data Preparation
  # ===========================================================================
  if (verbose) message("Preparing data...")

  keep_vars <- unique(c(
    id_col, time_col, y_col, r_col,
    covars_visit_XV, covars_outcome_fixed_XY, covars_outcome_random_link_ZY,
    covars_obs_fixed_XO, if (!is.null(covars_obs_random_link_ZO)) covars_obs_random_link_ZO
  ))
  df <- data[, keep_vars, drop = FALSE]

  n_row    <- nrow(df)
  id_vec   <- df[[id_col]]
  time_vec <- df[[time_col]]
  y_vec    <- df[[y_col]]
  r_vec    <- df[[r_col]]

  visit_mask <- !is.na(time_vec)

  subj_first_idx <- !duplicated(id_vec)
  id_levels      <- id_vec[subj_first_idx]
  n_subj         <- length(id_levels)
  subj_index     <- match(id_vec, id_levels)

  mi_visits_per_subj <- tabulate(subj_index[visit_mask], nbins = n_subj)

  # ===========================================================================
  # Stage 1: Visiting Parameter Estimation
  # ===========================================================================
  if (verbose) message("Stage 1: Estimating visiting process parameters...")

  XV_subject <- as.matrix(df[subj_first_idx, covars_visit_XV, drop = FALSE])

  visit_rows_idx <- which(visit_mask)
  n_visits_total <- length(visit_rows_idx)

  XV_visits_only <- as.matrix(df[visit_rows_idx, covars_visit_XV, drop = FALSE])
  V_target_sum   <- colSums(XV_visits_only, na.rm = TRUE)

  # Newton-Raphson for gamma
  gamma <- rep(0, ncol(XV_subject))
  for (iter in 1:100) {
    risk_score <- as.vector(exp(XV_subject %*% gamma))
    sum_risk <- sum(risk_score)

    XV_bar <- colSums(risk_score * XV_subject) / sum_risk
    score <- V_target_sum - XV_bar * n_visits_total

    XV_centered <- sweep(XV_subject, 2, XV_bar, "-")
    H <- -n_visits_total * crossprod(XV_centered * sqrt(risk_score)) / sum_risk

    step <- tryCatch(solve(H, score), error = function(e) rep(0, length(gamma)))
    gamma_new <- gamma - step

    if (max(abs(step)) < 1e-8) break
    gamma <- gamma_new
  }
  gamma_V_hat <- gamma
  names(gamma_V_hat) <- colnames(XV_subject)

  risk_scores_final <- as.vector(exp(XV_subject %*% gamma_V_hat))
  denom_risk_set    <- sum(risk_scores_final)

  unique_visit_times <- sort(unique(time_vec[visit_rows_idx]))
  dN_counts <- tabulate(match(time_vec[visit_rows_idx], unique_visit_times),
                        nbins = length(unique_visit_times))

  dLambda0_hat <- dN_counts / denom_risk_set
  names(dLambda0_hat) <- unique_visit_times

  C_global <- max(time_vec[visit_rows_idx], na.rm = TRUE)
  Lambda0_at_C_scalar <- sum(dLambda0_hat[unique_visit_times <= C_global])
  Lambda0_at_C_subj   <- rep(Lambda0_at_C_scalar, n_subj)

  nu_i <- risk_scores_final * Lambda0_at_C_subj

  mm_numerator <- sum(mi_visits_per_subj^2 - mi_visits_per_subj - nu_i^2, na.rm = TRUE)
  mm_denominator <- sum(nu_i^2, na.rm = TRUE)

  var_eta_mom <- if (mm_denominator > 0) max(mm_numerator / mm_denominator, 0) else 0

  sigma2_zeta_hat <- if (var_eta_mom > 0) log1p(var_eta_mom) else 0
  mu0_zeta        <- -0.5 * sigma2_zeta_hat

  # Laplace Approximation - vectorized
  if (sigma2_zeta_hat > 0) {
    nu_safe <- pmax(nu_i, 1e-12)
    log_eta_mode <- log((mi_visits_per_subj + 0.5) / nu_safe)

    for (k in 1:30) {
      exp_term <- exp(log(nu_safe) + log_eta_mode)
      g <- mi_visits_per_subj - exp_term - (log_eta_mode - mu0_zeta) / sigma2_zeta_hat
      H <- -exp_term - 1 / sigma2_zeta_hat
      step <- g / H
      step <- pmax(pmin(step, 2), -2)
      log_eta_mode <- log_eta_mode - step
      if (max(abs(step)) < 1e-8) break
    }
  } else {
    log_eta_mode <- rep(mu0_zeta, n_subj)
  }

  nu_exp_mode <- exp(pmax(log(nu_i), -700) + log_eta_mode)
  var_log_eta_post <- if (sigma2_zeta_hat > 0) 1 / (nu_exp_mode + 1 / sigma2_zeta_hat) else rep(0, n_subj)

  mu_Ui_hat     <- if (sigma2_zeta_hat > 0) (log_eta_mode - mu0_zeta) / sqrt(sigma2_zeta_hat) else rep(0, n_subj)
  sigma2_Ui_hat <- if (sigma2_zeta_hat > 0) var_log_eta_post / sigma2_zeta_hat else rep(0, n_subj)

  mu_Ui_row     <- mu_Ui_hat[subj_index]
  sigma2_Ui_row <- sigma2_Ui_hat[subj_index]

  if (verbose) {
    message(sprintf("  - Estimated sigma2_zeta = %.4f", sigma2_zeta_hat))
    message(sprintf("  - %d subjects with %d total visits", n_subj, n_visits_total))
  }

  # ===========================================================================
  # Stage 2: Observation Parameter Estimation via OU + PCL
  # ===========================================================================
  if (verbose) message("Stage 2: Estimating observation process via PCL...")

  obs_rows_idx <- which(!is.na(r_vec) & visit_mask)

  form_XO <- stats::as.formula(paste("~", paste(covars_obs_fixed_XO, collapse = "+")))
  XO_all  <- stats::model.matrix(form_XO, data = df)
  XO_obs  <- XO_all[obs_rows_idx, , drop = FALSE]
  R_obs   <- r_vec[obs_rows_idx]

  has_random_slopes <- !is.null(covars_obs_random_link_ZO) && length(covars_obs_random_link_ZO) > 0
  is_frailty_zero   <- (var_eta_mom < 1e-6)

  fit_obs_random <- (has_random_slopes && !is_frailty_zero)

  if (fit_obs_random) {
    form_ZO <- stats::as.formula(paste("~", paste(covars_obs_random_link_ZO, collapse = "+")))
    ZO_all  <- stats::model.matrix(form_ZO, data = df)
    ZO_obs  <- ZO_all[obs_rows_idx, , drop = FALSE]
  } else {
    ZO_all <- NULL
  }

  if (!fit_obs_random) {
    fit_glm <- stats::glm(stats::reformulate(covars_obs_fixed_XO, response = r_col),
                          data = df[obs_rows_idx, , drop = FALSE],
                          family = stats::binomial(link = "probit"))
    coef_temp <- stats::coef(fit_glm)
    alpha_O_hat <- numeric(ncol(XO_all)); names(alpha_O_hat) <- colnames(XO_all)
    alpha_O_hat[names(coef_temp)] <- coef_temp

    delta_O_hat <- setNames(numeric(0), character(0))
    Sigma_q_hat <- setNames(numeric(0), character(0))
    phi_hat     <- NA_real_
    pcl_value   <- NA_real_
    pcl_convergence <- NA_integer_

  } else {
    # Build pair data structure
    obs_df <- data.frame(
      obs_idx = seq_along(obs_rows_idx),
      row_idx = obs_rows_idx,
      subj = subj_index[obs_rows_idx],
      time = time_vec[obs_rows_idx],
      R = R_obs
    )
    obs_df <- obs_df[order(obs_df$subj, obs_df$time), ]

    subj_counts <- tabulate(obs_df$subj, nbins = n_subj)
    subj_with_pairs <- which(subj_counts >= 2)

    if (length(subj_with_pairs) == 0) {
      stop("No subjects with >= 2 visits for PCL estimation.", call. = FALSE)
    }

    # Pre-allocate pair vectors
    if (pair_type == "adjacent") {
      n_pairs_total <- sum(pmax(subj_counts - 1, 0))
    } else {
      n_pairs_total <- sum(choose(pmax(subj_counts, 0), 2))
    }

    idx_j <- integer(n_pairs_total)
    idx_k <- integer(n_pairs_total)

    ptr <- 1L
    subj_start <- c(0L, cumsum(subj_counts))

    for (s in subj_with_pairs) {
      start_idx <- subj_start[s] + 1L
      end_idx <- subj_start[s + 1]
      n_vis <- end_idx - start_idx + 1L

      if (n_vis < 2) next

      obs_indices <- start_idx:end_idx

      if (pair_type == "adjacent") {
        n_new_pairs <- n_vis - 1L
        idx_j[ptr:(ptr + n_new_pairs - 1L)] <- obs_indices[1:(n_vis - 1)]
        idx_k[ptr:(ptr + n_new_pairs - 1L)] <- obs_indices[2:n_vis]
      } else {
        pairs <- combn(n_vis, 2)
        n_new_pairs <- ncol(pairs)
        idx_j[ptr:(ptr + n_new_pairs - 1L)] <- obs_indices[pairs[1, ]]
        idx_k[ptr:(ptr + n_new_pairs - 1L)] <- obs_indices[pairs[2, ]]
      }
      ptr <- ptr + n_new_pairs
    }

    idx_j <- idx_j[1:(ptr - 1)]
    idx_k <- idx_k[1:(ptr - 1)]
    n_pairs <- length(idx_j)

    if (verbose) message(sprintf("  - Using %d %s pairs from %d subjects",
                                 n_pairs, pair_type, length(subj_with_pairs)))

    # Pre-extract ALL data needed for PCL (stored in parent environment)
    XO_j <- XO_obs[obs_df$obs_idx[idx_j], , drop = FALSE]
    XO_k <- XO_obs[obs_df$obs_idx[idx_k], , drop = FALSE]
    ZO_j <- ZO_obs[obs_df$obs_idx[idx_j], , drop = FALSE]
    ZO_k <- ZO_obs[obs_df$obs_idx[idx_k], , drop = FALSE]

    row_j <- obs_df$row_idx[idx_j]
    row_k <- obs_df$row_idx[idx_k]

    mu_U_pair <- mu_Ui_row[row_j]
    s2_U_pair <- sigma2_Ui_row[row_j]

    R_j <- obs_df$R[idx_j]
    R_k <- obs_df$R[idx_k]
    D_j <- 2 * R_j - 1
    D_k <- 2 * R_k - 1
    D_jk <- D_j * D_k

    delta_t <- abs(obs_df$time[idx_k] - obs_df$time[idx_j])

    # Pre-compute element-wise products for efficiency
    ZO_j_sq <- ZO_j^2
    ZO_k_sq <- ZO_k^2
    ZO_jk <- ZO_j * ZO_k

    p_alpha <- ncol(XO_obs)
    q_dim <- ncol(ZO_obs)

    # =========================================================================
    # PCL Negative Log-Likelihood (optimized)
    # =========================================================================
    neg_log_lik_PCL <- function(par) {
      alpha <- par[1:p_alpha]
      delta <- par[(p_alpha + 1):(p_alpha + q_dim)]
      sigma2_q <- exp(2 * par[(p_alpha + q_dim + 1):(p_alpha + 2 * q_dim)])
      phi <- exp(par[p_alpha + 2 * q_dim + 1])

      # Linear predictors
      lin_fixed_j <- as.vector(XO_j %*% alpha)
      lin_fixed_k <- as.vector(XO_k %*% alpha)
      delta_Z_j <- as.vector(ZO_j %*% delta)
      delta_Z_k <- as.vector(ZO_k %*% delta)

      M_j <- lin_fixed_j + delta_Z_j * mu_U_pair
      M_k <- lin_fixed_k + delta_Z_k * mu_U_pair

      # Variances
      ZSZ_j <- as.vector(ZO_j_sq %*% sigma2_q)
      ZSZ_k <- as.vector(ZO_k_sq %*% sigma2_q)

      V_j <- (delta_Z_j^2) * s2_U_pair + ZSZ_j + 1
      V_k <- (delta_Z_k^2) * s2_U_pair + ZSZ_k + 1

      # Covariance
      ZjSZk <- as.vector(ZO_jk %*% sigma2_q)
      rho_raw <- exp(-phi * delta_t)
      C_jk <- (delta_Z_j * delta_Z_k) * s2_U_pair + ZjSZk + rho_raw

      # Correlation
      inv_sqrt_Vj <- 1 / sqrt(pmax(V_j, 1e-10))
      inv_sqrt_Vk <- 1 / sqrt(pmax(V_k, 1e-10))
      Psi_jk <- pmax(pmin(C_jk * inv_sqrt_Vj * inv_sqrt_Vk, 0.9999), -0.9999)

      # Standardized arguments
      z_j <- pmax(pmin(D_j * M_j * inv_sqrt_Vj, 8), -8)
      z_k <- pmax(pmin(D_k * M_k * inv_sqrt_Vk, 8), -8)
      rho_jk <- D_jk * Psi_jk

      # Bivariate normal CDF
      prob <- pbvnorm_fast(z_j, z_k, rho_jk)

      -sum(log(pmax(prob, 1e-300)))
    }

    # =========================================================================
    # Analytic Gradient of PCL (KEY OPTIMIZATION)
    # =========================================================================
    grad_neg_log_lik_PCL <- function(par) {
      alpha <- par[1:p_alpha]
      delta <- par[(p_alpha + 1):(p_alpha + q_dim)]
      log_sigma_q <- par[(p_alpha + q_dim + 1):(p_alpha + 2 * q_dim)]
      log_phi <- par[p_alpha + 2 * q_dim + 1]

      sigma2_q <- exp(2 * log_sigma_q)
      phi <- exp(log_phi)

      # ----- Forward pass (same as neg_log_lik) -----
      lin_fixed_j <- as.vector(XO_j %*% alpha)
      lin_fixed_k <- as.vector(XO_k %*% alpha)
      delta_Z_j <- as.vector(ZO_j %*% delta)
      delta_Z_k <- as.vector(ZO_k %*% delta)

      M_j <- lin_fixed_j + delta_Z_j * mu_U_pair
      M_k <- lin_fixed_k + delta_Z_k * mu_U_pair

      ZSZ_j <- as.vector(ZO_j_sq %*% sigma2_q)
      ZSZ_k <- as.vector(ZO_k_sq %*% sigma2_q)

      V_j <- (delta_Z_j^2) * s2_U_pair + ZSZ_j + 1
      V_k <- (delta_Z_k^2) * s2_U_pair + ZSZ_k + 1

      ZjSZk <- as.vector(ZO_jk %*% sigma2_q)
      rho_raw <- exp(-phi * delta_t)
      C_jk <- (delta_Z_j * delta_Z_k) * s2_U_pair + ZjSZk + rho_raw

      inv_sqrt_Vj <- 1 / sqrt(pmax(V_j, 1e-10))
      inv_sqrt_Vk <- 1 / sqrt(pmax(V_k, 1e-10))
      Psi_jk <- C_jk * inv_sqrt_Vj * inv_sqrt_Vk
      Psi_jk_clip <- pmax(pmin(Psi_jk, 0.9999), -0.9999)

      z_j <- D_j * M_j * inv_sqrt_Vj
      z_k <- D_k * M_k * inv_sqrt_Vk
      z_j_clip <- pmax(pmin(z_j, 8), -8)
      z_k_clip <- pmax(pmin(z_k, 8), -8)
      rho_jk <- D_jk * Psi_jk_clip

      # Bivariate normal CDF
      prob <- pmax(pbvnorm_fast(z_j_clip, z_k_clip, rho_jk), 1e-300)

      # ----- Gradients of bivariate normal CDF -----
      # d/dz1 Phi2(z1,z2,rho) = phi(z1) * Phi((z2 - rho*z1)/sqrt(1-rho^2))
      # d/dz2 Phi2(z1,z2,rho) = phi(z2) * Phi((z1 - rho*z2)/sqrt(1-rho^2))
      # d/drho Phi2(z1,z2,rho) = phi2(z1,z2,rho) (bivariate density)

      one_minus_rho2 <- pmax(1 - rho_jk^2, 1e-10)
      sqrt_1_rho2 <- sqrt(one_minus_rho2)

      phi_zj <- dnorm(z_j_clip)
      phi_zk <- dnorm(z_k_clip)

      cond_arg_k <- (z_k_clip - rho_jk * z_j_clip) / sqrt_1_rho2
      cond_arg_j <- (z_j_clip - rho_jk * z_k_clip) / sqrt_1_rho2

      dprob_dzj <- phi_zj * pnorm(cond_arg_k)
      dprob_dzk <- phi_zk * pnorm(cond_arg_j)

      # Bivariate normal density
      quad_form <- (z_j_clip^2 + z_k_clip^2 - 2 * rho_jk * z_j_clip * z_k_clip) / one_minus_rho2
      dprob_drho <- exp(-0.5 * quad_form) / (2 * pi * sqrt_1_rho2)

      # ----- Chain rule: d(-log L)/d(param) = -sum(1/prob * dprob/dparam) -----
      w_j <- -dprob_dzj / prob
      w_k <- -dprob_dzk / prob
      w_rho <- -dprob_drho / prob

      # Derivatives of z_j, z_k, rho_jk w.r.t. intermediate variables
      # z_j = D_j * M_j / sqrt(V_j)
      # dz_j/dM_j = D_j / sqrt(V_j)
      # dz_j/dV_j = -D_j * M_j / (2 * V_j^{3/2}) = -z_j / (2*V_j)

      dL_dMj <- w_j * D_j * inv_sqrt_Vj
      dL_dMk <- w_k * D_k * inv_sqrt_Vk

      # rho_jk = D_jk * Psi_jk, Psi_jk = C_jk / (sqrt(V_j) * sqrt(V_k))
      # dPsi/dC_jk = 1 / (sqrt(V_j) * sqrt(V_k))
      # dPsi/dV_j = -Psi_jk / (2*V_j)
      # dPsi/dV_k = -Psi_jk / (2*V_k)

      w_rho_D <- w_rho * D_jk  # d(-logL)/dPsi_jk

      dL_dVj <- w_j * (-z_j / (2 * V_j)) + w_rho_D * (-Psi_jk / (2 * V_j))
      dL_dVk <- w_k * (-z_k / (2 * V_k)) + w_rho_D * (-Psi_jk / (2 * V_k))
      dL_dCjk <- w_rho_D * inv_sqrt_Vj * inv_sqrt_Vk

      # ----- Gradients w.r.t. parameters -----

      # alpha: dM_j/dalpha = XO_j, dM_k/dalpha = XO_k
      grad_alpha <- colSums(dL_dMj * XO_j) + colSums(dL_dMk * XO_k)

      # delta:
      # dM_j/ddelta = ZO_j * mu_U_pair
      # dV_j/ddelta = 2 * delta_Z_j * s2_U_pair * ZO_j
      # dC_jk/ddelta = (delta_Z_k * s2_U_pair) * ZO_j + (delta_Z_j * s2_U_pair) * ZO_k
      grad_delta <- colSums(dL_dMj * mu_U_pair * ZO_j) +
        colSums(dL_dMk * mu_U_pair * ZO_k) +
        colSums(dL_dVj * 2 * delta_Z_j * s2_U_pair * ZO_j) +
        colSums(dL_dVk * 2 * delta_Z_k * s2_U_pair * ZO_k) +
        colSums(dL_dCjk * delta_Z_k * s2_U_pair * ZO_j) +
        colSums(dL_dCjk * delta_Z_j * s2_U_pair * ZO_k)

      # log_sigma_q (chain rule: d(sigma2_q)/d(log_sigma_q) = 2*sigma2_q)
      # dV_j/d(sigma2_q[l]) = ZO_j[,l]^2
      # dV_k/d(sigma2_q[l]) = ZO_k[,l]^2
      # dC_jk/d(sigma2_q[l]) = ZO_j[,l] * ZO_k[,l]
      grad_log_sigma_q <- colSums(dL_dVj * ZO_j_sq) * 2 * sigma2_q +
        colSums(dL_dVk * ZO_k_sq) * 2 * sigma2_q +
        colSums(dL_dCjk * ZO_jk) * 2 * sigma2_q

      # log_phi: dC_jk/dphi = -delta_t * exp(-phi*delta_t), chain: dphi/d(log_phi) = phi
      grad_log_phi <- sum(dL_dCjk * (-delta_t * rho_raw)) * phi

      c(grad_alpha, grad_delta, grad_log_sigma_q, grad_log_phi)
    }

    # Initialize
    fit0 <- stats::glm(stats::reformulate(covars_obs_fixed_XO, response = r_col),
                       data = df[obs_rows_idx, , drop = FALSE],
                       family = stats::binomial(link = "probit"))
    alpha_init <- numeric(ncol(XO_obs)); names(alpha_init) <- colnames(XO_obs)
    tmp <- stats::coef(fit0); alpha_init[names(tmp)] <- tmp

    delta_init <- rep(0, q_dim)
    log_sigma_q_init <- rep(log(0.3), q_dim)
    log_phi_init <- log(phi_init)

    par_init <- c(alpha_init, delta_init, log_sigma_q_init, log_phi_init)

    lower_bounds <- c(rep(-Inf, p_alpha + q_dim), rep(log(1e-6), q_dim), log(phi_bounds[1]))
    upper_bounds <- c(rep(Inf, p_alpha + q_dim), rep(log(10), q_dim), log(phi_bounds[2]))

    # Optimize with analytic gradient
    opt <- stats::optim(
      par = par_init,
      fn = neg_log_lik_PCL,
      gr = grad_neg_log_lik_PCL,
      method = optim_method,
      lower = lower_bounds,
      upper = upper_bounds,
      control = list(maxit = 500, trace = if(verbose) 1 else 0, REPORT = 50)
    )

    if (opt$convergence != 0 && verbose) {
      warning(sprintf("PCL optimization may not have converged (code %d)", opt$convergence),
              call. = FALSE)
    }

    alpha_O_hat <- opt$par[1:p_alpha]
    names(alpha_O_hat) <- colnames(XO_obs)

    delta_O_hat <- opt$par[(p_alpha + 1):(p_alpha + q_dim)]
    names(delta_O_hat) <- colnames(ZO_obs)

    Sigma_q_hat <- exp(2 * opt$par[(p_alpha + q_dim + 1):(p_alpha + 2 * q_dim)])
    names(Sigma_q_hat) <- colnames(ZO_obs)

    phi_hat <- exp(opt$par[p_alpha + 2 * q_dim + 1])
    pcl_value <- opt$value
    pcl_convergence <- opt$convergence

    if (verbose) {
      # Compute median gap
      gaps <- numeric(0)
      for (s in seq_len(n_subj)) {
        idx <- which(subj_index == s & visit_mask)
        if (length(idx) >= 2) {
          times_s <- sort(time_vec[idx])
          gaps <- c(gaps, diff(times_s))
        }
      }
      median_gap <- if (length(gaps) > 0) median(gaps) else NA
      rho_median <- if (!is.na(median_gap)) exp(-phi_hat * median_gap) else NA
      message(sprintf("  - Estimated phi = %.4f (correlation at median gap %.2f: %.4f)",
                      phi_hat, median_gap, rho_median))
    }
  }

  # ===========================================================================
  # Calculate Marginal Probability omega_i(t) for ALL rows
  # ===========================================================================

  term_fixed_all <- as.vector(XO_all %*% alpha_O_hat)

  if (fit_obs_random) {
    term_latent_link_all <- as.vector(ZO_all %*% delta_O_hat)
    var_random_slope_all <- as.vector((ZO_all^2) %*% Sigma_q_hat)
  } else {
    term_latent_link_all <- numeric(n_row)
    var_random_slope_all <- numeric(n_row)
  }

  denom_scale_all <- sqrt(pmax(1 + var_random_slope_all + (term_latent_link_all^2) * sigma2_Ui_row, 1 + 1e-12))

  k_it_arg <- (term_fixed_all + term_latent_link_all * mu_Ui_row) / denom_scale_all
  k_it_clip <- pmax(pmin(k_it_arg, 12), -12)

  varpi_it_marginal_prob <- pnorm(k_it_clip)
  phi_k_it_density       <- dnorm(k_it_clip)

  # ===========================================================================
  # Stage 3: Outcome Parameter Estimation
  # ===========================================================================
  if (verbose) message("Stage 3: Estimating outcome model parameters...")

  correction_term <- (term_latent_link_all * sigma2_Ui_row / denom_scale_all) *
    (phi_k_it_density / pmax(varpi_it_marginal_prob, 1e-12))

  kappa_it <- mu_Ui_row + correction_term

  use_time_in_ZY <- time_col %in% covars_outcome_random_link_ZY
  ZY_vars_no_time <- setdiff(covars_outcome_random_link_ZY, time_col)

  if (length(ZY_vars_no_time) > 0) {
    ZY_matrix_row <- cbind(Intercept = 1, as.matrix(df[, ZY_vars_no_time, drop = FALSE]))
  } else {
    ZY_matrix_row <- matrix(1, nrow = n_row, ncol = 1, dimnames = list(NULL, "Intercept"))
  }
  if (use_time_in_ZY) ZY_matrix_row <- cbind(ZY_matrix_row, time = df[[time_col]])

  B_it_compensation <- sweep(ZY_matrix_row, 1, kappa_it, "*")

  # Weights
  XO_subj <- XO_all[subj_first_idx, , drop = FALSE]
  term_fixed_subj <- as.vector(XO_subj %*% alpha_O_hat)

  if (fit_obs_random) {
    ZO_subj <- ZO_all[subj_first_idx, , drop = FALSE]
    term_latent_link_subj <- as.vector(ZO_subj %*% delta_O_hat)
    var_random_slope_subj <- as.vector((ZO_subj^2) %*% Sigma_q_hat)
  } else {
    term_latent_link_subj <- numeric(n_subj)
    var_random_slope_subj <- numeric(n_subj)
  }

  denom_scale_subj <- sqrt(pmax(1 + var_random_slope_subj + (term_latent_link_subj^2) * sigma2_Ui_hat, 1 + 1e-12))
  k_it_subj <- (term_fixed_subj + term_latent_link_subj * mu_Ui_hat) / denom_scale_subj

  varpi_subj_base <- pnorm(pmax(pmin(k_it_subj, 12), -12))
  phi_subj_base   <- dnorm(pmax(pmin(k_it_subj, 12), -12))

  corr_subj <- (term_latent_link_subj * sigma2_Ui_hat / denom_scale_subj) *
    (phi_subj_base / pmax(varpi_subj_base, 1e-12))
  kappa_subj_base <- mu_Ui_hat + corr_subj

  weights_subj <- (mi_visits_per_subj * varpi_subj_base) / pmax(Lambda0_at_C_subj, 1e-12)
  sum_weights  <- sum(weights_subj)

  # Risk-Set Centering
  form_XY <- stats::as.formula(paste("~ 0 +", paste(covars_outcome_fixed_XY, collapse = "+")))
  XY_all  <- stats::model.matrix(form_XY, data = df)
  XY_subj <- XY_all[subj_first_idx, , drop = FALSE]

  XY_bar_vector <- colSums(XY_subj * weights_subj) / sum_weights
  XY_bar_matrix <- matrix(rep(XY_bar_vector, each = n_row), nrow = n_row,
                          dimnames = list(NULL, colnames(XY_subj)))

  if (length(ZY_vars_no_time) > 0) {
    ZY_subj <- cbind(Intercept = 1, as.matrix(df[subj_first_idx, ZY_vars_no_time, drop = FALSE]))
  } else {
    ZY_subj <- matrix(1, nrow = n_subj, ncol = 1, dimnames = list(NULL, "Intercept"))
  }
  B_subj <- sweep(ZY_subj, 1, kappa_subj_base, "*")

  B_bar_vector <- colSums(B_subj * weights_subj) / sum_weights
  B_bar_matrix <- matrix(rep(B_bar_vector, each = n_row), nrow = n_row,
                         dimnames = list(NULL, colnames(ZY_subj)))

  if (use_time_in_ZY) {
    kappa_bar <- sum(weights_subj * kappa_subj_base) / sum_weights
    B_bar_matrix <- cbind(B_bar_matrix, time = kappa_bar * df[[time_col]])
  }

  # WLS
  rows_event_outcome <- which(visit_mask & (r_vec == 1) & !is.na(y_vec))
  if (length(rows_event_outcome) == 0L) {
    stop("No observed outcome rows.", call. = FALSE)
  }

  Design_Matrix     <- cbind(XY_all, B_it_compensation)
  Design_Matrix_Bar <- cbind(XY_bar_matrix, B_bar_matrix)

  H_sub <- (Design_Matrix - Design_Matrix_Bar)[rows_event_outcome, , drop = FALSE]
  D_sub <- Design_Matrix[rows_event_outcome, , drop = FALSE]
  Y_sub <- y_vec[rows_event_outcome]

  # Use crossprod for efficiency
  LHS <- crossprod(H_sub, D_sub)  # = t(H_sub) %*% D_sub
  RHS <- crossprod(H_sub, Y_sub)  # = t(H_sub) %*% Y_sub

  coefs_final <- tryCatch(solve(LHS, RHS), error = function(e) base::qr.solve(LHS, RHS))
  coefs_final <- as.vector(coefs_final)

  p_xy <- ncol(XY_all)
  beta_Y_hat  <- coefs_final[1:p_xy]
  theta_Y_hat <- coefs_final[-(1:p_xy)]

  names(beta_Y_hat)  <- colnames(XY_all)
  names(theta_Y_hat) <- colnames(B_it_compensation)

  # Variance Estimation (Sandwich)
  res <- Y_sub - as.vector(D_sub %*% coefs_final)

  Bread <- crossprod(H_sub, D_sub)  # = t(H_sub) %*% D_sub
  score_contrib_i <- rowsum(H_sub * res, group = subj_index[rows_event_outcome])
  Meat <- crossprod(score_contrib_i)  # = t(S) %*% S where S is stacked scores

  invBread <- tryCatch(solve(Bread), error = function(e) base::qr.solve(Bread))

  # Sandwich: (H'D)^{-1} * Meat * (H'D)^{-T}
  V_robust <- invBread %*% Meat %*% t(invBread)

  se_beta <- sqrt(pmax(diag(V_robust)[1:p_xy], 0))
  names(se_beta) <- names(beta_Y_hat)

  if (verbose) {
    message(sprintf("  - Estimated %d fixed effects and %d latent link parameters",
                    length(beta_Y_hat), length(theta_Y_hat)))
  }

  # ===========================================================================
  # Compute median gap for output
  # ===========================================================================
  gaps <- numeric(0)
  for (s in seq_len(n_subj)) {
    idx <- which(subj_index == s & visit_mask)
    if (length(idx) >= 2) {
      times_s <- sort(time_vec[idx])
      gaps <- c(gaps, diff(times_s))
    }
  }
  median_gap <- if (length(gaps) > 0) median(gaps) else NA_real_

  # ===========================================================================
  # Return Results
  # ===========================================================================
  if (verbose) message("Done.")

  obs_params <- list(alpha_fixed = alpha_O_hat)
  if (fit_obs_random) {
    obs_params$delta_latent_random_link <- delta_O_hat
    obs_params$Sigma_q           <- Sigma_q_hat
    obs_params$phi_OU            <- phi_hat
    obs_params$pcl_value         <- pcl_value
    obs_params$convergence       <- pcl_convergence
  }

  structure(
    list(
      visiting_process = list(
        gamma_V           = gamma_V_hat,
        Lambda0_hat       = dLambda0_hat,
        sigma2_zeta       = sigma2_zeta_hat,
        mu_Ui             = setNames(mu_Ui_hat, id_levels),
        sigma2_Ui         = setNames(sigma2_Ui_hat, id_levels)
      ),
      observation_process = obs_params,
      longitudinal_outcome = list(
        beta_Y            = beta_Y_hat,
        theta_Y           = theta_Y_hat,
        se_beta           = se_beta,
        vcov_robust       = V_robust
      ),
      model_info = list(
        n_subjects  = n_subj,
        n_visits    = n_visits_total,
        n_observed  = sum(r_vec == 1, na.rm = TRUE),
        pair_type   = pair_type,
        median_gap  = median_gap
      ),
      call = call
    ),
    class = "CIMEHR_timevarying_ou"
  )
}

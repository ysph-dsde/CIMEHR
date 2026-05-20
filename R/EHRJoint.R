#' Joint Estimation for Electronic Health Record (EHR) Longitudinal
#' Processes (EHRJoint)
#'
#' @description
#' Fits the three-stage joint estimator of Du, Shi, and Mukherjee (2025) for
#' irregularly observed longitudinal Electronic Health Record (EHR) data with
#' informative visiting and informative observation processes. Stage 1
#' estimates the visiting-process intensity by solving an Andersen-Gill
#' partial-likelihood score equation. Stage 2 fits the observation indicator
#' via a probit (default) or logistic Generalized Linear Model (GLM). Stage 3
#' estimates the longitudinal outcome coefficients from a Liang-type
#' risk-set-centred estimating equation, augmented with a frailty-compensation
#' covariate constructed from the Stage 1 estimates. Analytic standard errors
#' are provided for all three stages: model-based for the visiting-process
#' partial likelihood and for the observation-process GLM, and a clustered
#' sandwich estimator for the outcome stage.
#'
#' @param data A \code{data.frame} in long (panel) format with one row per
#'   subject-time combination. Must contain columns for subject identifier,
#'   visit time, outcome, observation indicator, and all specified
#'   covariates. Rows with \code{NA} in the time column are treated as
#'   placeholder rows for zero-visit subjects.
#' @param id Character string giving the name of the subject identifier
#'   column. Default is \code{"id"}.
#' @param time Character string giving the name of the visit-time column.
#'   Default is \code{"time"}.
#' @param outcome Character string giving the name of the outcome column.
#'   Default is \code{"Y"}.
#' @param obs_indicator Character string giving the name of the binary
#'   observation indicator column (\code{1} = outcome observed, \code{0} =
#'   outcome missing at visit). Default is \code{"R"}.
#' @param visit_covs Character vector of covariate column names entering both
#'   the visiting-process intensity and the observation-process GLM. Default
#'   is \code{c("Z", "X")}.
#' @param long_covs Character vector of covariate column names for the
#'   fixed-effect part of the longitudinal outcome model. Default is
#'   \code{c("Z", "X")}.
#' @param random_covs Character vector of covariate column names that
#'   interact with the latent factor in the outcome model (frailty
#'   compensation). An intercept column is prepended automatically. Default
#'   is \code{"Z"}.
#' @param obs_link Character string giving the link function for the
#'   observation-process GLM. One of \code{"probit"} (default, matching the
#'   common data-generating process) or \code{"logit"}.
#' @param sigma_sq_floor Positive numeric scalar; numerical threshold below
#'   which the estimated frailty variance is treated as zero. When the
#'   estimate falls below this floor, the frailty-compensation covariates
#'   are dropped and the outcome model reduces to the no-correction
#'   submodel. Default is \code{1e-6}.
#'
#' @return An object of class \code{"EHRJoint"}, a named list with three
#'   components, each itself a named list:
#'
#' \describe{
#'   \item{\code{visiting_process}}{A list with elements:
#'     \describe{
#'       \item{\code{gamma_hat}}{Named numeric vector of visiting-process
#'         regression coefficients.}
#'       \item{\code{se_gamma}}{Named numeric vector of model-based standard
#'         errors for \code{gamma_hat}, from the inverse observed
#'         information matrix.}
#'       \item{\code{vcov_gamma}}{Numeric matrix; the model-based
#'         variance-covariance matrix for \code{gamma_hat}.}
#'       \item{\code{converged}}{Logical; \code{TRUE} if
#'         \code{\link[nleqslv]{nleqslv}} returned a successful termination
#'         code.}
#'     }
#'   }
#'   \item{\code{observation_process}}{A list with elements:
#'     \describe{
#'       \item{\code{alpha_hat}}{Named numeric vector of observation-process
#'         GLM coefficients (including intercept).}
#'       \item{\code{se_alpha}}{Named numeric vector of model-based standard
#'         errors for \code{alpha_hat}.}
#'       \item{\code{vcov_alpha}}{Numeric matrix; the model-based
#'         variance-covariance matrix for \code{alpha_hat}.}
#'       \item{\code{link}}{Character string; the link function used
#'         (\code{"probit"} or \code{"logit"}).}
#'     }
#'   }
#'   \item{\code{longitudinal_outcome}}{A list with elements:
#'     \describe{
#'       \item{\code{beta_hat}}{Named numeric vector of outcome fixed-effect
#'         coefficients.}
#'       \item{\code{theta_hat}}{Named numeric vector of latent-link
#'         coefficients (the part of the augmented coefficient vector
#'         multiplying the frailty-compensation covariates). All \code{NA}
#'         when the frailty variance is below \code{sigma_sq_floor}.}
#'       \item{\code{se_beta}}{Named numeric vector of clustered sandwich
#'         standard errors for \code{beta_hat}.}
#'       \item{\code{se_theta}}{Named numeric vector of clustered sandwich
#'         standard errors for \code{theta_hat}. All \code{NA} when the
#'         frailty variance is below \code{sigma_sq_floor}.}
#'       \item{\code{vcov_joint}}{Numeric matrix; the full clustered
#'         sandwich variance-covariance matrix for the augmented coefficient
#'         vector \code{c(beta_hat, theta_hat)}.}
#'       \item{\code{Lambda_C}}{Numeric scalar; Breslow estimate of the
#'         cumulative baseline visiting hazard over the full follow-up
#'         window.}
#'       \item{\code{sigma_sq}}{Numeric scalar; method-of-moments estimate
#'         of the frailty variance \eqn{\sigma^2}, truncated below at zero.}
#'       \item{\code{sigma_sq_raw}}{Numeric scalar; the untruncated
#'         method-of-moments estimate of the frailty variance.}
#'       \item{\code{use_B}}{Logical; \code{TRUE} if the frailty-compensation
#'         covariates were retained in the outcome model
#'         (\code{sigma_sq > sigma_sq_floor}), \code{FALSE} otherwise.}
#'       \item{\code{cond_LHS}}{Numeric scalar; the L2 condition number of
#'         the Stage 3 estimating-equation left-hand-side matrix. Values
#'         larger than \code{1e10} trigger a warning and a pseudoinverse
#'         fallback.}
#'     }
#'   }
#'   \item{\code{call}}{The matched call.}
#' }
#'
#' @details
#' The estimator proceeds in three stages.
#'
#' \strong{Stage 1 - Visiting Process.} The visit intensity is
#' \deqn{\lambda_i(t) = \eta_i \exp(\gamma^\top V_i)\,\lambda_0(t),}
#' where \eqn{\eta_i} is a shared frailty. The coefficient vector
#' \eqn{\gamma} is estimated by solving the Andersen-Gill partial-likelihood
#' score equation
#' \deqn{\sum_i \int [V_i - \bar V(\gamma; t)]\,dN_i(t) = 0,}
#' where \eqn{\bar V(\gamma; t)} is the risk-set-weighted mean of \eqn{V}.
#' Only visit rows enter the score; placeholder rows for zero-visit subjects
#' are excluded. Standard errors are obtained from the inverse of the
#' observed information matrix.
#'
#' \strong{Stage 2 - Observation Process.} The observation indicator is
#' modelled by
#' \deqn{P(R_{ij} = 1 \mid V_i) = g^{-1}(\alpha^\top V_i^*),}
#' where \eqn{V_i^* = (1, V_i)} and \eqn{g} is the probit (default) or
#' logistic link. The coefficient vector \eqn{\alpha} is estimated by
#' \code{\link[stats]{glm}} and the model-based variance-covariance matrix
#' is returned.
#'
#' \strong{Stage 3 - Longitudinal Outcome.} An empirical Bayes estimator of
#' the latent loading \eqn{B_i} is constructed as
#' \deqn{\hat B_i = \frac{1 + n_i \hat\sigma^2}{1 +
#'   \hat\Lambda_i \hat\sigma^2} - 1,}
#' where \eqn{n_i} is the visit count and
#' \eqn{\hat\Lambda_i = \exp(\hat\gamma^\top V_i) \hat\Lambda_0(C)}.
#' The outcome model
#' \deqn{Y_i(t) = \beta^\top X_i(t) + \theta^\top \hat B_i Q_i^*(t) +
#'   \varepsilon_i(t)}
#' is estimated by solving the risk-set-centred estimating equation
#' \deqn{\sum_{i,j} (M_{ij} - \bar M)(Y_{ij} - M_{ij}^\top \beta^*) = 0,}
#' with \eqn{M_{ij} = (X_{ij}, \hat B_i Q_{ij}^*)} and the centering
#' weighted by \eqn{n_i \hat w_i} where \eqn{\hat w_i} is the predicted
#' observation probability. The equation is linear in \eqn{\beta^*} and is
#' solved directly by \code{\link[base]{solve}}, with a Moore-Penrose
#' pseudoinverse fallback (via \code{\link[MASS]{ginv}}) on singular systems.
#' A clustered sandwich estimator
#' \deqn{\hat V = B^{-1} \Bigl(\sum_i A_i A_i^\top\Bigr) B^{-\top},
#'   \quad B = H^\top M,\
#'   A_i = \sum_j H_{ij}(Y_{ij} - M_{ij}^\top \hat\beta^*)}
#' is returned. This variance is conditional on the Stage 1 nuisance
#' estimates and does not propagate uncertainty from
#' \eqn{(\hat\gamma, \hat\sigma^2, \hat\Lambda_0)}.
#'
#' If the estimated frailty variance \eqn{\hat\sigma^2} falls below
#' \code{sigma_sq_floor}, the frailty-compensation block is dropped and the
#' outcome model reduces to the no-correction submodel with only
#' \code{long_covs}; \code{theta_hat} and \code{se_theta} are then returned
#' as \code{NA}.
#'
#' @note
#' This function requires the \pkg{nleqslv} and \pkg{MASS} packages.
#'
#' @references
#' Du, J., Shi, X., and Mukherjee, B. (2025). A new statistical approach for
#' joint modeling of longitudinal outcomes measured in electronic health
#' records with clinically informative presence and observation processes.
#' arXiv:2410.13113. \doi{10.48550/arXiv.2410.13113}.
#'
#' @seealso
#' \code{\link{Joint_modeling_visiting_and_longitudinal_Liang}} for a related
#' joint model without an observation-process component;
#' \code{\link{CIMEHR}} for the three-stage CIMEHR estimator with a full
#' shared latent structure;
#' \code{\link{bootstrap}} for bootstrap supplementation of the analytic
#' standard errors.
#'
#' @examples
#' dat <- sim_ehr_data[sim_ehr_data$id <= 30, ]
#' fit <- EHRJoint(
#'   dat,
#'   visit_covs  = "Age",
#'   long_covs   = "Age",
#'   random_covs = "Age"
#' )
#' print(fit)
#'
#' @importFrom stats as.formula binomial coef glm plogis pnorm vcov
#' @export
EHRJoint <- function(data,
                     id = "id",
                     time = "time",
                     outcome = "Y",
                     obs_indicator = "R",
                     visit_covs = c("Z", "X"),
                     long_covs = c("Z", "X"),
                     random_covs = "Z",
                     obs_link = c("probit", "logit"),
                     sigma_sq_floor = 1e-6) {

  if (!requireNamespace("nleqslv", quietly = TRUE)) {
    stop("Package 'nleqslv' is required. Install it with: install.packages('nleqslv')",
         call. = FALSE)
  }
  if (!requireNamespace("MASS", quietly = TRUE)) {
    stop("Package 'MASS' is required. Install it with: install.packages('MASS')",
         call. = FALSE)
  }
  obs_link <- match.arg(obs_link)

  call <- match.call()

  visit_covs  <- .check_covars(visit_covs,  data, "visit_covs")
  long_covs   <- .check_covars(long_covs,   data, "long_covs")
  random_covs <- .check_covars(random_covs, data, "random_covs")
  if (length(visit_covs)  == 0L) stop("`visit_covs` must include at least one covariate.",  call. = FALSE)
  if (length(long_covs)   == 0L) stop("`long_covs` must include at least one covariate.",   call. = FALSE)
  if (length(random_covs) == 0L) stop("`random_covs` must include at least one covariate.", call. = FALSE)

  # ------------------------------------------------------------------
  # Data prep
  # ------------------------------------------------------------------
  id_vec     <- data[[id]]
  unique_ids <- unique(id_vec)
  n_subj     <- length(unique_ids)
  subj_index <- match(id_vec, unique_ids)
  sf         <- which(!duplicated(id_vec))            # baseline-row indices

  V  <- as.matrix(data[, visit_covs, drop = FALSE])
  VV <- V[sf, , drop = FALSE]                         # one row per subject

  X  <- as.matrix(data[, long_covs, drop = FALSE])
  XX <- X[sf, , drop = FALSE]

  Q_full  <- as.matrix(data[, random_covs, drop = FALSE])
  Q0      <- cbind(`(Intercept)` = 1, Q_full)         # nrow(data) x (1 + p_Q)
  Q0_base <- Q0[sf, , drop = FALSE]

  Time <- data[[time]]
  y    <- data[[outcome]]

  visit_mask     <- !is.na(Time)
  n_visits_total <- sum(visit_mask)
  ni <- tabulate(subj_index[visit_mask], nbins = n_subj)

  # ------------------------------------------------------------------
  # Stage 1: Andersen-Gill (AG) partial-likelihood score for gamma
  # (sum only over visit rows; exclude placeholder rows for zero-visit
  # subjects)
  # ------------------------------------------------------------------
  V_visits     <- V[visit_mask, , drop = FALSE]
  sum_V_visits <- colSums(V_visits)                   # sum_i N_i V_i

  gamma_score <- function(gamma) {
    eta   <- as.vector(exp(VV %*% gamma))
    denom <- sum(eta)
    numer <- colSums(VV * eta)
    sum_V_visits - (numer / denom) * n_visits_total
  }

  g0 <- rep(0, ncol(V))
  gamma_solve <- nleqslv::nleqslv(g0, gamma_score, control = list(trace = 0))
  if (!gamma_solve$termcd %in% c(1L, 2L)) {
    warning("Stage 1 (gamma) did not converge: termcd = ", gamma_solve$termcd)
  }
  gamma_hat <- gamma_solve$x
  names(gamma_hat) <- visit_covs

  exp_gW   <- as.vector(exp(VV %*% gamma_hat))        # length n_subj
  V_denom  <- sum(exp_gW)
  V_bar    <- colSums(VV * exp_gW) / V_denom
  VV_cent  <- sweep(VV, 2, V_bar, "-")
  Info_g   <- n_visits_total * crossprod(VV_cent * sqrt(exp_gW)) / V_denom
  vcov_g   <- tryCatch(solve(Info_g),
                       error = function(e) MASS::ginv(Info_g))
  se_gamma <- sqrt(pmax(diag(vcov_g), 0))
  names(se_gamma) <- visit_covs
  dimnames(vcov_g) <- list(visit_covs, visit_covs)

  # ------------------------------------------------------------------
  # Stage 2: observation-process GLM (probit by default)
  # ------------------------------------------------------------------
  obs_formula <- stats::as.formula(
    paste(obs_indicator, "~", paste(visit_covs, collapse = " + "))
  )
  fit_IP <- stats::glm(obs_formula, data = data,
                       family = stats::binomial(link = obs_link))
  alpha_hat  <- stats::coef(fit_IP)
  vcov_alpha <- stats::vcov(fit_IP)
  se_alpha   <- sqrt(pmax(diag(vcov_alpha), 0))
  names(se_alpha) <- names(alpha_hat)

  link_inv <- if (obs_link == "logit") stats::plogis else stats::pnorm
  V_int    <- cbind(`(Intercept)` = 1, V)
  wi_full  <- as.vector(link_inv(V_int %*% alpha_hat))
  wi_base  <- wi_full[sf]                             # length n_subj

  # ------------------------------------------------------------------
  # Stage 3a: Lambda_C (Breslow) and sigma^2 (method of moments)
  # ------------------------------------------------------------------
  Lambda_C_est <- n_visits_total / V_denom
  Lambda_i     <- exp_gW * Lambda_C_est               # length n_subj
  sigma_sq_raw <- sum(ni^2 - ni - Lambda_i^2) / sum(Lambda_i^2)
  sigma_sq_hat <- max(sigma_sq_raw, 0)

  # ------------------------------------------------------------------
  # Stage 3b: frailty-compensation B_hat (broadcast by subject index)
  # ------------------------------------------------------------------
  use_B <- (sigma_sq_hat > sigma_sq_floor)
  if (use_B) {
    B_i       <- (1 + ni * sigma_sq_hat) / (1 + Lambda_i * sigma_sq_hat) - 1
    Bhat      <- Q0      * B_i[subj_index]            # broadcast by subject
    Bhat_base <- Q0_base * B_i                        # n_subj x (1 + p_Q)
  } else {
    B_i       <- rep(0, n_subj)
    Bhat      <- matrix(0, nrow = nrow(data), ncol = ncol(Q0),
                        dimnames = list(NULL, colnames(Q0)))
    Bhat_base <- matrix(0, nrow = n_subj,    ncol = ncol(Q0),
                        dimnames = list(NULL, colnames(Q0)))
  }

  # ------------------------------------------------------------------
  # Stage 3c: Liang-style risk-set-centred estimating equation,
  # solved as a direct linear system
  # ------------------------------------------------------------------
  w_subj  <- ni * wi_base                             # E[M_i] proxy
  denom_w <- sum(w_subj)
  if (denom_w <= 0) stop("Stage 3: total subject weight is zero.", call. = FALSE)

  Xbar <- as.vector(colSums(w_subj * XX) / denom_w)
  if (use_B) {
    Bbar      <- as.vector(colSums(w_subj * Bhat_base) / denom_w)
    design_M  <- cbind(X, Bhat)
    centering <- c(Xbar, Bbar)
  } else {
    design_M  <- X
    centering <- Xbar
  }

  Dbar <- matrix(centering, nrow = nrow(data), ncol = length(centering),
                 byrow = TRUE)
  H    <- design_M - Dbar

  obs_mask <- !is.na(y)
  H_obs <- H[obs_mask, , drop = FALSE]
  D_obs <- design_M[obs_mask, , drop = FALSE]
  Y_obs <- y[obs_mask]

  LHS <- crossprod(H_obs, D_obs)
  RHS <- as.vector(crossprod(H_obs, Y_obs))
  cond_LHS <- kappa(LHS)
  if (cond_LHS > 1e10) {
    warning("Stage 3: LHS condition number = ", signif(cond_LHS, 3),
            "; estimates may be unstable.", call. = FALSE)
  }
  coefs_all <- tryCatch(
    as.vector(solve(LHS, RHS)),
    error = function(e) as.vector(MASS::ginv(LHS) %*% RHS)
  )

  p_X <- ncol(X)
  p_B <- if (use_B) ncol(Bhat) else 0
  beta_hat  <- coefs_all[seq_len(p_X)]
  theta_hat <- if (use_B) coefs_all[p_X + seq_len(p_B)] else
    rep(NA_real_, ncol(Q0))
  names(beta_hat)  <- long_covs
  names(theta_hat) <- colnames(Q0)

  # Subject-clustered sandwich variance
  res_obs  <- Y_obs - as.vector(D_obs %*% coefs_all)
  score_ij <- H_obs * res_obs
  score_i  <- rowsum(score_ij, subj_index[obs_mask])
  Meat     <- crossprod(score_i)
  invA     <- tryCatch(solve(LHS), error = function(e) MASS::ginv(LHS))
  Vsand    <- invA %*% Meat %*% t(invA)

  se_all   <- sqrt(pmax(diag(Vsand), 0))
  se_beta  <- se_all[seq_len(p_X)]
  se_theta <- if (use_B) se_all[p_X + seq_len(p_B)] else
    rep(NA_real_, ncol(Q0))
  names(se_beta)  <- long_covs
  names(se_theta) <- colnames(Q0)

  vcov_names <- c(long_covs, if (use_B) colnames(Q0))
  if (length(vcov_names) == ncol(Vsand))
    dimnames(Vsand) <- list(vcov_names, vcov_names)

  structure(list(
    visiting_process = list(
      gamma_hat  = gamma_hat,
      se_gamma   = se_gamma,
      vcov_gamma = vcov_g,
      converged  = (gamma_solve$termcd %in% c(1L, 2L))
    ),
    observation_process = list(
      alpha_hat  = alpha_hat,
      se_alpha   = se_alpha,
      vcov_alpha = vcov_alpha,
      link       = obs_link
    ),
    longitudinal_outcome = list(
      beta_hat     = beta_hat,
      theta_hat    = theta_hat,
      se_beta      = se_beta,
      se_theta     = se_theta,
      vcov_joint   = Vsand,
      Lambda_C     = Lambda_C_est,
      sigma_sq     = sigma_sq_hat,
      sigma_sq_raw = sigma_sq_raw,
      use_B        = use_B,
      cond_LHS     = cond_LHS
    ),
    call = call
  ), class = "EHRJoint")
}

#' Generalized Informative Visiting and Examination in EHR (Gauss-Hermite Version)
#'
#' Fits a joint model for visiting, observation, and longitudinal processes using
#' Gauss-Hermite quadrature and online filtering for time-varying observation probabilities.
#'
#' @param data Data frame containing longitudinal data.
#' @param id_col Subject ID column name. Default is "id".
#' @param time_col Time variable column name. Default is "time".
#' @param y_col Outcome variable column name. Default is "Y".
#' @param r_col Observation indicator column name. Default is "R".
#' @param r_lag_col Lagged observation indicator column name. If NULL, no time-varying lambda.
#' @param covars_visit_XV Covariates for visiting process. Default is c("Z", "X").
#' @param covars_outcome_fixed_XY Fixed effect covariates for outcome. Default is c("Z", "X").
#' @param covars_outcome_link_ZY Link covariates for outcome (random effects). Default is c("Z").
#' @param covars_obs_fixed_XO Fixed effect covariates for observation process. Default is c("Z", "X").
#' @param covars_obs_random_ZO Random effect covariates for observation process. Default is "Z".
#' @param gh_points_U Number of Gauss-Hermite points for U integration. Default is 7.
#' @param gh_points_q Number of Gauss-Hermite points for q integration. Default is 5.
#' @param gh_clip_eta Clipping value for eta in probit. Default is 12.
#' @param clip_k Clipping value for k in Phi computation. Default is 12.
#' @param filter_max_iter Maximum iterations for Laplace filter. Default is 6.
#' @param verbose If TRUE, print progress messages. Default is FALSE.
#'
#' @return A list containing estimates from visiting, observation, and longitudinal processes.
#'
#' @export
GIVEHR_GH <- function(
    data,
    id_col                  = "id",
    time_col                = "time",
    y_col                   = "Y",
    r_col                   = "R",
    r_lag_col               = NULL,
    covars_visit_XV         = c("Z", "X"),
    covars_outcome_fixed_XY = c("Z", "X"),
    covars_outcome_link_ZY  = c("Z"),
    covars_obs_fixed_XO     = c("Z", "X"),
    covars_obs_random_ZO    = "Z",
    gh_points_U             = 7,
    gh_points_q             = 5,
    gh_clip_eta             = 12,
    clip_k                  = 12,
    filter_max_iter         = 6,
    verbose                 = FALSE
) {

  if (!requireNamespace("nleqslv", quietly = TRUE)) stop("Package 'nleqslv' is required.")

  # ==========================================================================
  # Helper Functions
  # ==========================================================================

  # Gauss-Hermite quadrature nodes and weights
  .ghq <- function(n) {
    i <- 1:(n - 1)
    b <- sqrt(i / 2)
    J <- diag(0, n, n)
    J[cbind(1:(n - 1), 2:n)] <- b
    J[cbind(2:n, 1:(n - 1))] <- b
    eg <- eigen(J, symmetric = TRUE)
    ord <- order(eg$values)
    list(x = eg$values[ord], w = (eg$vectors[1, ]^2)[ord] * sqrt(pi))
  }

  # Stable row-wise log-sum-exp
  row_logsumexp <- function(mat) {
    n <- nrow(mat)
    K <- ncol(mat)
    if (K == 1L) return(mat[, 1L])
    out <- numeric(n)
    for (i in seq_len(n)) {
      row_i <- mat[i, ]
      m <- max(row_i)
      if (!is.finite(m)) { out[i] <- m; next }
      out[i] <- m + log(sum(exp(row_i - m)))
    }
    out
  }

  log_add_exp_vec <- function(x, y) {
    m <- pmax(x, y)
    m + log(exp(x - m) + exp(y - m))
  }

  # Compute omega, kappa for a single subject given full posterior (m, V)
  compute_omega_kappa_single <- function(m_vec, V_mat, a_const, aw, clip_k) {
    mu_L <- a_const + sum(aw * m_vec)
    var_L <- as.numeric(t(aw) %*% V_mat %*% aw)
    var_L <- max(var_L, 0)

    denom <- sqrt(1 + var_L)
    k <- mu_L / denom
    k <- max(min(k, clip_k), -clip_k)

    Phi_k <- pnorm(k)
    Phi_k <- min(max(Phi_k, 1e-12), 1 - 1e-12)
    phi_k <- dnorm(k)

    omega <- Phi_k
    cov_UL <- sum(V_mat[1, ] * aw)
    kappa <- m_vec[1] + (cov_UL / denom) * (phi_k / Phi_k)

    list(omega = omega, kappa = kappa)
  }

  # Joint Laplace update for w = (U, eps) with probit likelihood
  update_joint_laplace <- function(m_prior, V_prior, R, a_const, aw,
                                   max_iter, tol = 1e-8, clip_eta) {
    p_dim <- length(m_prior)

    P <- tryCatch(solve(V_prior), error = function(e) {
      diag(1 / pmax(diag(V_prior), 1e-12))
    })

    w <- m_prior

    for (it in 1:max_iter) {
      eta <- a_const + sum(aw * w)
      eta <- max(min(eta, clip_eta), -clip_eta)

      prob <- pnorm(eta)
      prob <- min(max(prob, 1e-12), 1 - 1e-12)
      phi <- dnorm(eta)

      if (R == 1) {
        imr <- phi / prob
        d1eta <- imr
        d2eta <- -imr * (eta + imr)
      } else {
        p0 <- 1 - prob
        p0 <- min(max(p0, 1e-12), 1 - 1e-12)
        imr_u <- phi / p0
        d1eta <- -imr_u
        d2eta <- -imr_u * (-eta + imr_u)
      }

      g <- -as.vector(P %*% (w - m_prior)) + aw * d1eta
      H <- -P + d2eta * tcrossprod(aw)
      H <- 0.5 * (H + t(H))

      step <- tryCatch(solve(H, g), error = function(e) qr.solve(H, g))
      step <- pmax(pmin(step, 1), -1)

      w_new <- w - step
      if (max(abs(w_new - w)) < tol) { w <- w_new; break }
      w <- w_new
    }

    # Final Hessian at mode
    eta <- a_const + sum(aw * w)
    eta <- max(min(eta, clip_eta), -clip_eta)
    prob <- pnorm(eta)
    prob <- min(max(prob, 1e-12), 1 - 1e-12)
    phi <- dnorm(eta)

    if (R == 1) {
      imr <- phi / prob
      d2eta <- -imr * (eta + imr)
    } else {
      p0 <- 1 - prob
      p0 <- min(max(p0, 1e-12), 1 - 1e-12)
      imr_u <- phi / p0
      d2eta <- -imr_u * (-eta + imr_u)
    }

    H <- -P + d2eta * tcrossprod(aw)
    H <- 0.5 * (H + t(H))
    V_post <- -tryCatch(solve(H), error = function(e) qr.solve(H))

    eig <- eigen(V_post, symmetric = TRUE)
    eig$values <- pmax(eig$values, 1e-10)
    V_post <- eig$vectors %*% diag(eig$values) %*% t(eig$vectors)

    list(m = w, V = V_post)
  }

  # ==========================================================================
  # 0. Data Preparation
  # ==========================================================================
  stopifnot(all(c(id_col, time_col, y_col, r_col) %in% names(data)))

  df <- data
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

  mi_visits_per_subj <- tabulate(subj_index[visit_mask], n_subj)

  # Time-varying omega flag

  has_time_varying_omega <- !is.null(r_lag_col)

  # Internal R_lag: first visit per subject = 0
  r_lag_vec <- rep(0L, n_row)
  if (has_time_varying_omega) {
    v_idx <- which(visit_mask)
    o <- order(subj_index[v_idx], time_vec[v_idx])
    v_sorted <- v_idx[o]
    sid_sorted <- subj_index[v_sorted]
    r_sorted <- r_vec[v_sorted]

    lag_sorted <- c(0L, head(r_sorted, -1))
    new_subj <- c(TRUE, sid_sorted[-1] != sid_sorted[-length(sid_sorted)])
    lag_sorted[new_subj] <- 0L
    r_lag_vec[v_sorted] <- lag_sorted
  }

  # ==========================================================================
  # Stage 1: Visiting Process
  # ==========================================================================
  visit_rows_idx <- which(visit_mask)
  XV_subject <- as.matrix(df[subj_first_idx, covars_visit_XV, drop = FALSE])
  XV_visits  <- as.matrix(df[visit_rows_idx, covars_visit_XV, drop = FALSE])

  score_gamma <- function(g) {
    w <- exp(XV_subject %*% g)
    sw <- sum(w)
    colSums(XV_visits) - as.vector(crossprod(w, XV_subject)) / sw * length(visit_rows_idx)
  }

  gamma_V_hat <- nleqslv::nleqslv(rep(0, ncol(XV_subject)), score_gamma)$x
  names(gamma_V_hat) <- colnames(XV_subject)

  exp_gW <- as.vector(exp(XV_subject %*% gamma_V_hat))
  Lambda0_C <- length(visit_rows_idx) / sum(exp_gW)

  nu_i <- exp_gW * Lambda0_C
  var_eta_mom <- max((sum(mi_visits_per_subj^2 - mi_visits_per_subj) - sum(nu_i^2)) / sum(nu_i^2), 0)
  sigma2_zeta_hat <- if (var_eta_mom > 0) log1p(var_eta_mom) else 0
  mu0_zeta <- -0.5 * sigma2_zeta_hat
  sigma_zeta <- sqrt(sigma2_zeta_hat)

  # EB posterior for U
  if (sigma2_zeta_hat > 0) {
    log_eta_mode <- log(pmax(mi_visits_per_subj, 0.5) / pmax(nu_i, 1e-12))
    for (iter in 1:20) {
      em <- exp(log_eta_mode)
      g <- mi_visits_per_subj - nu_i * em - (log_eta_mode - mu0_zeta) / sigma2_zeta_hat
      H <- -nu_i * em - 1 / sigma2_zeta_hat
      log_eta_mode <- log_eta_mode - pmin(pmax(g / H, -1), 1)
    }
    log_eta_var <- 1 / (nu_i * exp(log_eta_mode) + 1 / sigma2_zeta_hat)
  } else {
    log_eta_mode <- rep(mu0_zeta, n_subj)
    log_eta_var <- rep(0, n_subj)
  }

  mu_Ui_hat <- if (sigma2_zeta_hat > 0) (log_eta_mode - mu0_zeta) / sigma_zeta else rep(0, n_subj)
  sigma2_Ui_hat <- if (sigma2_zeta_hat > 0) log_eta_var / sigma2_zeta_hat else rep(0, n_subj)

  # ==========================================================================
  # Stage 2: Observation Process (GH Quadrature)
  # ==========================================================================
  form_XO <- as.formula(paste("~", paste(covars_obs_fixed_XO, collapse = "+")))
  XO_all <- model.matrix(form_XO, data = df)

  d <- 0
  ZO_all <- NULL
  if (!is.null(covars_obs_random_ZO) && length(covars_obs_random_ZO) > 0) {
    form_ZO <- as.formula(paste("~", paste(covars_obs_random_ZO, collapse = "+")))
    ZO_all <- model.matrix(form_ZO, data = df)
    d <- ncol(ZO_all)
  }

  obs_rows_idx <- which(visit_mask & !is.na(r_vec))
  ord <- order(subj_index[obs_rows_idx], time_vec[obs_rows_idx])
  obs_ord <- obs_rows_idx[ord]

  R_obs <- r_vec[obs_ord]
  Rlag_obs <- if (has_time_varying_omega) r_lag_vec[obs_ord] else rep(0, length(obs_ord))
  XO_obs <- XO_all[obs_ord, , drop = FALSE]
  ZO_obs <- if (d > 0) ZO_all[obs_ord, , drop = FALSE] else matrix(0, length(obs_ord), 0)
  sid_obs <- subj_index[obs_ord]

  rle_sid <- rle(sid_obs)
  sid_unique <- rle_sid$values
  subj_m <- rle_sid$lengths
  n_obs_subj <- length(sid_unique)
  obs_to_subj <- rep(seq_len(n_obs_subj), subj_m)

  p <- ncol(XO_all)

  # GH quadrature setup
  gh_U <- .ghq(gh_points_U)
  K_U <- length(gh_U$x)
  log_wU <- log(gh_U$w / sqrt(pi))

  if (d > 0) {
    gh_q <- .ghq(gh_points_q)
    K_q <- length(gh_q$x)
    log_wq1 <- log(gh_q$w / sqrt(pi))

    if (d == 1) {
      q_grid_base <- matrix(gh_q$x, ncol = 1)
      log_wq_grid <- log_wq1
    } else {
      idx_mat <- as.matrix(do.call(expand.grid, rep(list(seq_len(K_q)), d)))
      q_grid_base <- matrix(gh_q$x[idx_mat], ncol = d)
      log_wq_grid <- rowSums(matrix(log_wq1[idx_mat], ncol = d))
    }
    n_q <- nrow(q_grid_base)
  } else {
    q_grid_base <- matrix(0, 1, 0)
    log_wq_grid <- 0
    n_q <- 1
  }

  mu_subj <- mu_Ui_hat[sid_unique]
  sd_subj <- sqrt(2 * pmax(sigma2_Ui_hat[sid_unique], 0))
  U_nodes_subj <- mu_subj + outer(sd_subj, gh_U$x)

  is1 <- (R_obs == 1)

  # Negative log-likelihood for Stage 2
  negloglik <- function(par) {
    alpha <- par[1:p]

    if (d > 0) {
      delta <- par[(p + 1):(p + d)]
      log_sq <- par[(p + d + 1):(p + 2 * d)]
      sq <- exp(log_sq)

      Q <- sqrt(2) * sweep(q_grid_base, 2, sq, "*")
      S_q_contrib <- ZO_obs %*% t(Q)
      S_delta <- as.vector(ZO_obs %*% delta)
    } else {
      delta <- numeric(0)
      sq <- numeric(0)
      S_q_contrib <- matrix(0, nrow(XO_obs), 1)
      S_delta <- rep(0, nrow(XO_obs))
    }

    lambda <- if (has_time_varying_omega) par[length(par)] else 0
    Va <- as.vector(XO_obs %*% alpha)

    log_total <- rep(-Inf, n_obs_subj)

    for (k in 1:K_U) {
      U_obs_k <- U_nodes_subj[obs_to_subj, k]
      eta_base <- Va + S_delta * U_obs_k
      if (has_time_varying_omega) eta_base <- eta_base + lambda * Rlag_obs

      eta_mat <- sweep(S_q_contrib, 1, eta_base, "+")
      eta_mat <- pmax(pmin(eta_mat, gh_clip_eta), -gh_clip_eta)

      logPhi <- pnorm(eta_mat, log.p = TRUE)
      log1mPhi <- pnorm(eta_mat, lower.tail = FALSE, log.p = TRUE)

      log_p_mat <- log1mPhi
      log_p_mat[is1, ] <- logPhi[is1, , drop = FALSE]

      log_lik_subj_q <- rowsum(log_p_mat, obs_to_subj, reorder = FALSE)
      if (n_q > 1) {
        log_lik_subj_q <- sweep(log_lik_subj_q, 2, log_wq_grid, "+")
        log_q_sum <- row_logsumexp(log_lik_subj_q)
      } else {
        log_q_sum <- as.vector(log_lik_subj_q[, 1] + log_wq_grid[1])
      }

      log_term <- log_q_sum + log_wU[k]
      log_total <- log_add_exp_vec(log_total, log_term)
    }

    -sum(log_total)
  }

  # Initialize and optimize
  fit0 <- glm.fit(XO_obs, R_obs, family = binomial("probit"))
  alpha0 <- coef(fit0)
  alpha0[is.na(alpha0)] <- 0
  par0 <- alpha0
  if (d > 0) par0 <- c(par0, rep(0, d), rep(-2, d))
  if (has_time_varying_omega) par0 <- c(par0, 0)

  opt <- nlminb(par0, negloglik, control = list(eval.max = 350, iter.max = 250, rel.tol = 1e-6))
  if (verbose) message("Stage2: conv=", opt$convergence, " obj=", round(opt$objective, 3))

  alpha_O_hat <- opt$par[1:p]
  names(alpha_O_hat) <- colnames(XO_all)

  if (d > 0) {
    delta_O_hat <- opt$par[(p + 1):(p + d)]
    names(delta_O_hat) <- colnames(ZO_all)
    sigma_q_hat <- exp(opt$par[(p + d + 1):(p + 2 * d)])
    names(sigma_q_hat) <- colnames(ZO_all)
    Sigma_q_hat <- diag(sigma_q_hat^2, d, d)
  } else {
    delta_O_hat <- numeric(0)
    sigma_q_hat <- numeric(0)
    Sigma_q_hat <- matrix(0, 0, 0)
  }

  lambda_hat <- if (has_time_varying_omega) opt$par[length(opt$par)] else 0

  # ==========================================================================
  # Stage 3: Outcome (Joint Filtering with Online-Sum)
  # ==========================================================================

  # Outcome design matrices
  form_XY <- as.formula(paste("~ 0 +", paste(covars_outcome_fixed_XY, collapse = "+")))
  XY_all <- model.matrix(form_XY, data = df)
  XY_subject <- XY_all[subj_first_idx, , drop = FALSE]
  px <- ncol(XY_all)

  use_time_in_ZY <- time_col %in% covars_outcome_link_ZY
  ZY_vars_no_time <- setdiff(covars_outcome_link_ZY, time_col)

  df_subj <- df[subj_first_idx, , drop = FALSE]
  if (length(ZY_vars_no_time) > 0) {
    ZY_subject_base <- cbind(Intercept = 1, as.matrix(df_subj[, ZY_vars_no_time, drop = FALSE]))
    ZY_row_base <- cbind(Intercept = 1, as.matrix(df[, ZY_vars_no_time, drop = FALSE]))
  } else {
    ZY_subject_base <- matrix(1, n_subj, 1, dimnames = list(NULL, "Intercept"))
    ZY_row_base <- matrix(1, n_row, 1, dimnames = list(NULL, "Intercept"))
  }
  ZY_row <- if (use_time_in_ZY) cbind(ZY_row_base, time = df[[time_col]]) else ZY_row_base
  pz <- ncol(ZY_subject_base)

  # Observation model quantities at subject level
  XO_subject <- XO_all[subj_first_idx, , drop = FALSE]
  ZO_subject <- if (d > 0) ZO_all[subj_first_idx, , drop = FALSE] else matrix(0, n_subj, 0)
  a_subj_base <- as.vector(XO_subject %*% alpha_O_hat)
  b_subj <- if (d > 0) as.vector(ZO_subject %*% delta_O_hat) else rep(0, n_subj)

  # Weight factor
  weight_factor <- mi_visits_per_subj / pmax(Lambda0_C, 1e-12)

  # Dimension of joint latent vector w = (U, eps)
  dim_w <- 1 + d

  # Initialize posterior states
  m_post <- vector("list", n_subj)
  V_post <- vector("list", n_subj)

  for (i in seq_len(n_subj)) {
    m_post[[i]] <- c(mu_Ui_hat[i], rep(0, d))
    V_i <- matrix(0, dim_w, dim_w)
    V_i[1, 1] <- sigma2_Ui_hat[i]
    if (d > 0) diag(V_i)[2:(1 + d)] <- diag(Sigma_q_hat)
    V_post[[i]] <- V_i
  }

  R_prev_all <- rep(0L, n_subj)

  omega_cur <- numeric(n_subj)
  kappa_cur <- numeric(n_subj)
  p_cur <- numeric(n_subj)

  get_aw <- function(k) {
    if (d > 0) c(b_subj[k], as.numeric(ZO_subject[k, ])) else b_subj[k]
  }

  # Initialize omega/kappa
  for (i in seq_len(n_subj)) {
    a_const <- a_subj_base[i] + lambda_hat * R_prev_all[i]
    aw <- get_aw(i)
    ok <- compute_omega_kappa_single(m_post[[i]], V_post[[i]], a_const, aw, clip_k)
    omega_cur[i] <- ok$omega
    kappa_cur[i] <- ok$kappa
    p_cur[i] <- omega_cur[i] * weight_factor[i]
  }

  # Initialize global sums
  S_p <- sum(p_cur)
  S_pX <- colSums(p_cur * XY_subject)
  S_pkZ <- colSums(p_cur * kappa_cur * ZY_subject_base)
  S_pk <- sum(p_cur * kappa_cur)

  # Event times
  event_times <- sort(unique(time_vec[visit_mask]))
  nT <- length(event_times)

  # Build visit structure
  v_idx <- which(visit_mask)
  visit_time_idx <- match(time_vec[v_idx], event_times)
  visit_subj_idx <- subj_index[v_idx]
  visit_r <- r_vec[v_idx]

  time_order <- order(visit_time_idx, visit_subj_idx)
  v_idx_sorted <- v_idx[time_order]
  visit_time_sorted <- visit_time_idx[time_order]
  visit_subj_sorted <- visit_subj_idx[time_order]
  visit_r_sorted <- visit_r[time_order]

  time_rle <- rle(visit_time_sorted)
  time_unique <- time_rle$values
  time_lengths <- time_rle$lengths
  time_ends <- cumsum(time_lengths)
  time_starts <- c(1L, time_ends[-length(time_ends)] + 1L)
  names(time_starts) <- as.character(time_unique)
  names(time_ends) <- as.character(time_unique)

  # Storage for row-level results
  XY_bar_row <- matrix(0, n_row, px)
  colnames(XY_bar_row) <- colnames(XY_all)

  B_bar_row <- matrix(0, n_row, ncol(ZY_row))
  colnames(B_bar_row) <- colnames(ZY_row)

  kappa_row <- rep(NA_real_, n_row)

  # Process event times sequentially
  for (t_idx in seq_len(nT)) {
    XY_bar_t <- S_pX / max(S_p, 1e-12)
    B_bar_base_t <- S_pkZ / max(S_p, 1e-12)
    kappa_bar_t <- S_pk / max(S_p, 1e-12)

    t_char <- as.character(t_idx)
    if (t_char %in% names(time_starts)) {
      idx_start <- time_starts[t_char]
      idx_end <- time_ends[t_char]

      for (j in idx_start:idx_end) {
        row_idx <- v_idx_sorted[j]
        k <- visit_subj_sorted[j]
        R_obs_k <- visit_r_sorted[j]

        XY_bar_row[row_idx, ] <- XY_bar_t
        B_bar_row[row_idx, colnames(ZY_subject_base)] <- B_bar_base_t
        if (use_time_in_ZY) {
          B_bar_row[row_idx, "time"] <- kappa_bar_t * time_vec[row_idx]
        }

        kappa_row[row_idx] <- kappa_cur[k]

        if (is.na(R_obs_k)) next

        # Subtract old contributions
        S_p <- S_p - p_cur[k]
        S_pX <- S_pX - p_cur[k] * XY_subject[k, ]
        S_pkZ <- S_pkZ - p_cur[k] * kappa_cur[k] * ZY_subject_base[k, ]
        S_pk <- S_pk - p_cur[k] * kappa_cur[k]

        # Update posterior
        a_const <- a_subj_base[k] + lambda_hat * R_prev_all[k]
        aw <- get_aw(k)

        upd <- update_joint_laplace(
          m_prior = m_post[[k]], V_prior = V_post[[k]],
          R = R_obs_k, a_const = a_const, aw = aw,
          max_iter = filter_max_iter, clip_eta = gh_clip_eta
        )
        m_post[[k]] <- upd$m
        V_post[[k]] <- upd$V

        R_prev_all[k] <- R_obs_k

        # Recompute omega, kappa, p
        a_const_new <- a_subj_base[k] + lambda_hat * R_prev_all[k]
        ok <- compute_omega_kappa_single(m_post[[k]], V_post[[k]], a_const_new, aw, clip_k)
        omega_cur[k] <- ok$omega
        kappa_cur[k] <- ok$kappa
        p_cur[k] <- omega_cur[k] * weight_factor[k]

        # Add new contributions
        S_p <- S_p + p_cur[k]
        S_pX <- S_pX + p_cur[k] * XY_subject[k, ]
        S_pkZ <- S_pkZ + p_cur[k] * kappa_cur[k] * ZY_subject_base[k, ]
        S_pk <- S_pk + p_cur[k] * kappa_cur[k]
      }
    }
  }

  # Build B_it_compensation
  B_it_compensation <- matrix(0, n_row, ncol(ZY_row))
  colnames(B_it_compensation) <- colnames(ZY_row)
  valid_kappa <- !is.na(kappa_row)
  B_it_compensation[valid_kappa, ] <- ZY_row[valid_kappa, , drop = FALSE] * kappa_row[valid_kappa]

  # Solve WLS
  rows_event_outcome <- which(visit_mask & r_vec == 1 & !is.na(y_vec) & !is.na(kappa_row))
  if (length(rows_event_outcome) == 0L) stop("No rows with valid time, R==1, Y observed, and kappa computed.")

  Design_Matrix <- cbind(XY_all, B_it_compensation)
  Design_Matrix_Bar <- cbind(XY_bar_row, B_bar_row)

  H_sub <- Design_Matrix[rows_event_outcome, , drop = FALSE] - Design_Matrix_Bar[rows_event_outcome, , drop = FALSE]
  D_sub <- Design_Matrix[rows_event_outcome, , drop = FALSE]
  Y_sub <- y_vec[rows_event_outcome]

  LHS <- crossprod(H_sub, D_sub)
  RHS <- crossprod(H_sub, Y_sub)

  coefs_final <- tryCatch(solve(LHS, RHS), error = function(e) qr.solve(LHS, RHS))

  beta_Y_hat <- as.vector(coefs_final[1:px])
  theta_Y_hat <- as.vector(coefs_final[-(1:px)])
  names(beta_Y_hat) <- colnames(XY_all)
  names(theta_Y_hat) <- colnames(ZY_row)

  # ==========================================================================
  # Return Results
  # ==========================================================================
  obs_params <- list(alpha_fixed = alpha_O_hat)
  if (d > 0) {
    obs_params$delta_latent_link <- delta_O_hat
    obs_params$Sigma_q <- Sigma_q_hat
  }
  if (has_time_varying_omega) {
    obs_params$lambda <- lambda_hat
  }

  list(
    visiting_process = list(
      gamma_V = gamma_V_hat,
      Lambda0_C = Lambda0_C,
      sigma2_zeta = sigma2_zeta_hat,
      mu_Ui = setNames(mu_Ui_hat, id_levels),
      sigma2_Ui = setNames(sigma2_Ui_hat, id_levels)
    ),
    observation_process = obs_params,
    longitudinal_outcome = list(
      beta_Y = beta_Y_hat,
      theta_Y = theta_Y_hat
    )
  )
}

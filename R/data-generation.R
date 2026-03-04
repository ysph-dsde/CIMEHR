# Simulate longitudinal data for GIVEHR with OU+PCL observation process
# Matches the model specification in Web Appendix S9
sim_data_gen <- function(
    n,
    time.start = 0,
    time.end   = 60,
    seed       = 123,

    # ---- Baseline covariates ----
    pZ  = 0.5,
    muX = 0, sdX = 1,

    # ---- Visiting process ----
    # lambda_i(t) = eta_i * exp(gamma' X^V) * lambda_0(t)
    # eta_i = exp(mu0 + sigma_zeta * U_i), U_i ~ N(0,1)
    gamma = c(intercept = -2.2, Z = 1, X = 1),
    sigma_zeta = 1.0,  # frailty SD on log scale

    # ---- Observation process (Probit with OU) ----
    # R_ij = I(R_tilde_ij > 0)
    # R_tilde_ij = alpha' X^O + q_i' Z^O + xi_i(t)
    # q_i | U_i ~ N(delta * U_i, Sigma_q)
    # xi_i(t) ~ OU process with Var=1, Cov(s,t) = exp(-phi|t-s|)
    alpha = c(intercept = 1.5, Z = -0.5, X = -0.3),  # fixed effects
    delta = c(intercept = -0.3, Z = -0.5),          # latent link (delta' Z^O)
    sigma_q = c(intercept = 0.3, Z = 0.2),          # sqrt of diag(Sigma_q)
    phi = 0.1,  # OU decay rate

    # ---- Outcome process ----
    # Y_i(t) = beta0(t) + beta' X^Y + b_i' Z^Y + epsilon
    # b_i | U_i ~ N(theta * U_i, Sigma_b)
    beta    = c(intercept = -2, Z = -0.5, X = 0.5),
    beta_t  = 0.1,
    theta   = c(intercept = 0.5, Z = 0.3),  # outcome link
    Sigma_b = matrix(c(1, 0, 0, 2), 2, 2),  # residual RE covariance
    sigma_y = 1
) {

  set.seed(seed)

  # ============================================================
  # Helper: Generate OU process at irregular times
  # ============================================================
  generate_OU <- function(times, phi) {
    n_times <- length(times)
    if (n_times == 0) return(numeric(0))
    if (n_times == 1) return(rnorm(1, 0, 1))

    # Sort times (should already be sorted, but ensure)
    ord <- order(times)
    times_sorted <- times[ord]

    xi <- numeric(n_times)
    xi[1] <- rnorm(1, 0, 1)  # Stationary distribution

    for (j in 2:n_times) {
      dt <- times_sorted[j] - times_sorted[j-1]
      rho <- exp(-phi * dt)
      xi[j] <- rho * xi[j-1] + sqrt(1 - rho^2) * rnorm(1)
    }

    # Return in original order
    xi[order(ord)]
  }

  # ============================================================
  # Helper: Homogeneous Poisson arrivals
  # ============================================================
  r_poisson_times <- function(lambda, max_time) {
    if (!is.finite(lambda) || lambda <= 0) return(numeric(0))
    t <- 0
    times <- numeric(0)
    repeat {
      t <- t + rexp(1, rate = lambda)
      if (t > max_time) break
      times <- c(times, t)
    }
    times
  }

  # ============================================================
  # 1. Generate baseline covariates
  # ============================================================
  base_df <- data.frame(
    id = seq_len(n),
    Z  = rbinom(n, 1, pZ),
    X  = rnorm(n, muX, sdX)
  )

  # ============================================================
  # 2. Generate shared latent factor U_i ~ N(0, 1)
  # ============================================================
  U_i <- rnorm(n, 0, 1)
  base_df$U <- U_i

  # ============================================================
  # 3. Generate frailty eta_i for visiting process
  # ============================================================
  mu0 <- -0.5 * sigma_zeta^2  # Ensure E[eta] = 1
  log_eta <- mu0 + sigma_zeta * U_i
  eta_i <- exp(log_eta)
  base_df$eta <- eta_i

  # ============================================================
  # 4. Generate visit times (Poisson process)
  # ============================================================
  tau <- time.end - time.start

  # Visit intensity: lambda_i = eta_i * exp(gamma' X^V)
  XV <- cbind(1, base_df$Z, base_df$X)
  gamma_vec <- c(gamma["intercept"], gamma["Z"], gamma["X"])
  lambda_i <- eta_i * exp(as.vector(XV %*% gamma_vec))

  visit_times_list <- lapply(lambda_i, function(lam) {
    r_poisson_times(lam, tau) + time.start
  })

  # Count visits per subject
  n_visits <- sapply(visit_times_list, length)
  base_df$n_visits <- n_visits

  # ============================================================
  # 5. Generate observation process random effects q_i
  # ============================================================
  # q_i | U_i ~ N(delta * U_i, Sigma_q)
  # q_i = delta * U_i + b_tilde_i, where b_tilde_i ~ N(0, Sigma_q)

  q_dim <- length(delta)
  delta_vec <- c(delta["intercept"], delta["Z"])
  sigma_q_vec <- c(sigma_q["intercept"], sigma_q["Z"])
  Sigma_q_diag <- sigma_q_vec^2

  # Residual random effects for observation process
  b_tilde_obs <- matrix(rnorm(n * q_dim), nrow = n, ncol = q_dim)
  b_tilde_obs <- sweep(b_tilde_obs, 2, sigma_q_vec, "*")  # Scale by sigma_q

  # Full q_i = delta * U_i + b_tilde
  q_i <- outer(U_i, delta_vec) + b_tilde_obs
  colnames(q_i) <- c("q_intercept", "q_Z")
  base_df <- cbind(base_df, q_i)

  # ============================================================
  # 6. Generate outcome process random effects b_i
  # ============================================================
  # b_i | U_i ~ N(theta * U_i, Sigma_b)

  theta_vec <- c(theta["intercept"], theta["Z"])

  # Residual random effects for outcome
  b_resid <- MASS::mvrnorm(n, mu = c(0, 0), Sigma = Sigma_b)

  # Full b_i = theta * U_i + b_resid
  b_i <- outer(U_i, theta_vec) + b_resid
  colnames(b_i) <- c("b_intercept", "b_Z")
  base_df <- cbind(base_df, b_i)

  # ============================================================
  # 7. Expand to long format
  # ============================================================
  long_list <- vector("list", n)

  for (i in seq_len(n)) {
    times_i <- visit_times_list[[i]]
    m_i <- length(times_i)

    if (m_i == 0) {
      # Subject with no visits
      long_list[[i]] <- data.frame(
        id = i,
        Z = base_df$Z[i],
        X = base_df$X[i],
        time = NA_real_,
        R = NA_integer_,
        Y = NA_real_,
        U = base_df$U[i],
        eta = base_df$eta[i],
        q_intercept = base_df$q_intercept[i],
        q_Z = base_df$q_Z[i],
        b_intercept = base_df$b_intercept[i],
        b_Z = base_df$b_Z[i]
      )
    } else {
      # Generate OU process for this subject
      xi_i <- generate_OU(times_i, phi)

      # Observation process design matrix Z^O
      ZO_i <- cbind(1, rep(base_df$Z[i], m_i))  # [1, Z]

      # Fixed effects design matrix X^O
      XO_i <- cbind(1, rep(base_df$Z[i], m_i), rep(base_df$X[i], m_i))  # [1, Z, X]

      # Linear predictor for observation
      alpha_vec <- c(alpha["intercept"], alpha["Z"], alpha["X"])
      lin_fixed <- as.vector(XO_i %*% alpha_vec)

      # Random effects part: q_i' Z^O
      q_i_vec <- c(base_df$q_intercept[i], base_df$q_Z[i])
      lin_random <- as.vector(ZO_i %*% q_i_vec)

      # Full latent propensity: R_tilde = alpha' X^O + q_i' Z^O + xi(t)
      R_tilde <- lin_fixed + lin_random + xi_i

      # Observation indicator
      R_i <- as.integer(R_tilde > 0)

      # Outcome model
      beta_vec <- c(beta["intercept"], beta["Z"], beta["X"])
      XY_i <- cbind(1, rep(base_df$Z[i], m_i), rep(base_df$X[i], m_i))

      # Fixed effects + time effect
      mu_Y <- as.vector(XY_i %*% beta_vec) + beta_t * times_i

      # Random effects: b_i' Z^Y where Z^Y = [1, Z]
      ZY_i <- cbind(1, rep(base_df$Z[i], m_i))
      b_i_vec <- c(base_df$b_intercept[i], base_df$b_Z[i])
      mu_Y <- mu_Y + as.vector(ZY_i %*% b_i_vec)

      # Add noise
      Y_i <- mu_Y + rnorm(m_i, 0, sigma_y)

      # Set Y to NA when R = 0
      Y_i[R_i == 0] <- NA_real_

      long_list[[i]] <- data.frame(
        id = i,
        Z = base_df$Z[i],
        X = base_df$X[i],
        time = times_i,
        R = R_i,
        Y = Y_i,
        U = base_df$U[i],
        eta = base_df$eta[i],
        q_intercept = base_df$q_intercept[i],
        q_Z = base_df$q_Z[i],
        b_intercept = base_df$b_intercept[i],
        b_Z = base_df$b_Z[i]
      )
    }
  }

  long_df <- do.call(rbind, long_list)
  rownames(long_df) <- NULL

  # ============================================================
  # 8. Summary statistics
  # ============================================================
  prop_no_visit <- mean(n_visits == 0)
  obs_rate <- mean(long_df$R[!is.na(long_df$R)])

  cat(sprintf("Simulation summary (n=%d):\n", n))
  cat(sprintf("  Subjects with no visits: %.1f%%\n", 100 * prop_no_visit))
  cat(sprintf("  Total visits: %d\n", sum(n_visits)))
  cat(sprintf("  Median visits per subject: %.1f\n", median(n_visits)))
  cat(sprintf("  Observation rate: %.1f%%\n", 100 * obs_rate))
  cat(sprintf("  Total observed outcomes: %d\n", sum(long_df$R == 1, na.rm = TRUE)))

  # ============================================================
  # Return
  # ============================================================
  structure(
    list(
      long_data = long_df,
      subject_data = base_df,
      prop_no_visit = prop_no_visit,
      obs_rate = obs_rate,
      params = list(
        n = n, time.start = time.start, time.end = time.end, seed = seed,
        pZ = pZ, muX = muX, sdX = sdX,
        gamma = gamma, sigma_zeta = sigma_zeta,
        alpha = alpha, delta = delta, sigma_q = sigma_q, phi = phi,
        beta = beta, beta_t = beta_t, theta = theta, Sigma_b = Sigma_b, sigma_y = sigma_y
      )
    ),
    class = "sim_GIVEHR_OU_PCL"
  )
}

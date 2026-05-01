#' Simulate Longitudinal Electronic Health Record (EHR) Data with Informative
#' Visiting and Observation Processes
#'
#' @description
#' Generates synthetic longitudinal data from a joint model that includes an
#' informative visiting process (Poisson with log-normal frailty), an
#' informative observation process (probit with latent-factor-linked random
#' effects), and a linear longitudinal outcome model with
#' latent-factor-linked random effects. By default, no Ornstein-Uhlenbeck (OU)
#' serial-correlation component is added; an optional OU component is enabled
#' by setting \code{phi > 0}.
#'
#' Covariates are user-specified through the \code{covariates} argument as a
#' named list, and per-process coefficients are passed as named numeric vectors
#' keyed by those covariate names. This lets the same generator produce
#' arbitrary covariate sets (e.g., generic Z/X, or EHR-style Age/Gender/NSES)
#' without changes to the function. See the Examples section.
#'
#' The data-generating mechanism matches the model specification in
#' Yang, Shi, and Mukherjee (2026).
#'
#' @param n Integer; number of subjects to simulate.
#' @param time.start Numeric scalar; start of the observation window. Default
#'   is \code{0}.
#' @param time.end Numeric scalar; end of the observation window (the constant
#'   censoring time used for all subjects). Default is \code{60}.
#' @param seed Integer; random seed for reproducibility. Default is \code{123}.
#' @param covariates Named list specifying the baseline covariates to generate.
#'   Each element is itself a list with at least a \code{type} entry, which
#'   must be either \code{"continuous"} or \code{"binary"}. Continuous
#'   covariates accept additional fields \code{dist} (\code{"normal"} or
#'   \code{"uniform"}; default \code{"normal"}), \code{mean}/\code{sd} (for
#'   normal), \code{min}/\code{max} (for uniform), and \code{standardize}
#'   (logical; whether to z-standardize after generation, default \code{FALSE}).
#'   Binary covariates accept \code{prob}, the Bernoulli success probability
#'   (default \code{0.5}). When \code{covariates = NULL} (the default), a
#'   minimal generic spec \code{list(Z = list(type = "binary", prob = 0.5),
#'   X = list(type = "continuous", dist = "normal", mean = 0, sd = 1))} is
#'   used so that the function works out of the box for the package's
#'   internal tests and existing examples.
#' @param gamma Named numeric vector of visiting-process log-intensity
#'   coefficients. Must include an \code{intercept} entry; remaining names
#'   should match a subset of the covariates declared in \code{covariates}.
#'   Covariates declared but not named in \code{gamma} are excluded from the
#'   visiting model (coefficient zero).
#' @param sigma_zeta Positive numeric scalar; standard deviation of the
#'   log-normal frailty on the log scale. Default is \code{1.0}.
#' @param alpha Named numeric vector of observation-process probit
#'   fixed-effect coefficients. Must include an \code{intercept}; remaining
#'   names should be a subset of \code{covariates}.
#' @param obs_random_covs Character vector naming the covariates whose
#'   observation-process random slopes are linked to the shared latent
#'   factor. An intercept random effect is always included and need not be
#'   listed here. The default is \code{"Z"}, matching the default
#'   \code{delta} and \code{sigma_q} dimensions; pass \code{character(0)}
#'   for an intercept-only random effect.
#' @param delta Named numeric vector of latent-link coefficients for the
#'   observation process. Must include an \code{intercept}; the remaining
#'   names should match \code{obs_random_covs}.
#' @param sigma_q Named numeric vector of residual standard deviations for
#'   the observation-process random effects. Must include an
#'   \code{intercept}; the remaining names should match
#'   \code{obs_random_covs}.
#' @param phi Nonnegative numeric scalar controlling the optional OU
#'   serial-correlation component in the observation process. Default
#'   \code{0} (no OU component).
#' @param beta Named numeric vector of outcome fixed-effect coefficients
#'   (excluding the time slope). Must include \code{intercept}; remaining
#'   names should be a subset of \code{covariates}.
#' @param beta_t Numeric scalar; coefficient for the linear time effect in
#'   the outcome model. Default is \code{0.1}.
#' @param outcome_random_covs Character vector naming the covariates whose
#'   outcome-model random slopes are linked to the shared latent factor. An
#'   intercept random effect is always included. The default is \code{"Z"}.
#' @param theta Named numeric vector of outcome latent-link coefficients.
#'   Must include \code{intercept}; remaining names should match
#'   \code{outcome_random_covs}.
#' @param Sigma_b Numeric positive-definite matrix; the residual covariance
#'   of the outcome random effects. Dimension must equal
#'   \code{1 + length(outcome_random_covs)}. Default is the
#'   \eqn{2 \times 2} matrix \code{matrix(c(1, 0, 0, 2), 2, 2)} for a single
#'   random slope.
#' @param sigma_y Positive numeric scalar; residual standard deviation of
#'   the outcome measurement error. Default is \code{1}.
#' @param verbose Logical; if \code{TRUE} (the default), print a short
#'   simulation summary via \code{\link{message}}.
#'
#' @return An object of class \code{"sim_CIMEHR"}, a named list with components
#'   \code{long_data} (the long-format \code{data.frame} returned to the user;
#'   contains \code{id}, the declared covariates, \code{time}, \code{R},
#'   \code{Y}, and \code{C}), \code{latent} (a subject-level
#'   \code{data.frame} with the latent factor \eqn{U_i}, frailty
#'   \eqn{\eta_i}, and random effects, retained for diagnostics and not
#'   merged into \code{long_data}), \code{prop_no_visit}, \code{obs_rate},
#'   and \code{params} (echoes all data-generating parameters for use as
#'   \code{true_values} in \code{method_comparisons()}).
#'
#' @section Latent quantities:
#' The latent variables \eqn{U_i}, \eqn{\eta_i}, \eqn{q_i}, and \eqn{b_i} are
#' not returned in \code{long_data} because they are unobservable in real EHR
#' data; including them would obscure the analyst's view of what a user
#' actually has. They remain available in the \code{latent} component for
#' simulation diagnostics.
#'
#' @seealso \code{\link{CIMEHR}}, \code{\link{CIMEHR_timevarying_integral}},
#'   \code{\link{CIMEHR_timevarying_ou}}.
#'
#' @examples
#' # Generic two-covariate setup (back-compatible with earlier package versions)
#' dat <- sim_data_gen(n = 30, seed = 42)
#' head(dat$long_data)
#'
#' # EHR-style multi-covariate setup
#' ehr_covs <- list(
#'   Age      = list(type = "continuous", dist = "uniform",
#'                   min = 18, max = 99, standardize = TRUE),
#'   Gender   = list(type = "binary", prob = 0.5),
#'   Marital  = list(type = "binary", prob = 0.4),
#'   Black    = list(type = "binary", prob = 0.30),
#'   Hispanic = list(type = "binary", prob = 0.20),
#'   NSES     = list(type = "continuous", dist = "normal",
#'                   mean = 0, sd = 1, standardize = FALSE)
#' )
#' dat_ehr <- sim_data_gen(
#'   n = 50, seed = 1, covariates = ehr_covs,
#'   gamma = c(intercept = -2.0, Age = 0.10, Marital = -0.05,
#'             Black = 0.15, Hispanic = -0.05, NSES = 0.12),
#'   alpha = c(intercept = 0.5, Age = 0.10, Gender = -0.10, NSES = 0.20),
#'   obs_random_covs     = "NSES",
#'   delta   = c(intercept = -0.3, NSES = -0.5),
#'   sigma_q = c(intercept =  0.3, NSES =  0.2),
#'   beta = c(intercept = -2.0, Age = 0.05, Gender = 0.10,
#'            Black = 0.15, Hispanic = 0.05, NSES = -0.25),
#'   outcome_random_covs = "NSES",
#'   theta   = c(intercept = 0.5, NSES = 0.3)
#' )
#' head(dat_ehr$long_data)
#'
#' @importFrom stats rbinom rnorm rexp runif median
#' @export
sim_data_gen <- function(
    n,
    time.start = 0,
    time.end   = 60,
    seed       = 123,

    # ---- Covariate specification ----
    covariates = NULL,

    # ---- Visiting process ----
    gamma      = c(intercept = -2.2, Z = 1, X = 1),
    sigma_zeta = 1.0,

    # ---- Observation process ----
    alpha           = c(intercept = 1.5, Z = -0.5, X = -0.3),
    obs_random_covs = "Z",
    delta           = c(intercept = -0.3, Z = -0.5),
    sigma_q         = c(intercept =  0.3, Z =  0.2),
    phi             = 0,

    # ---- Outcome process ----
    beta                = c(intercept = -2, Z = -0.5, X = 0.5),
    beta_t              = 0.1,
    outcome_random_covs = "Z",
    theta               = c(intercept = 0.5, Z = 0.3),
    Sigma_b             = matrix(c(1, 0, 0, 2), 2, 2),
    sigma_y             = 1,

    verbose = TRUE
) {

  # Save and restore RNG state so we don't mutate the user's global stream
  if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    on.exit(assign(".Random.seed", old_seed, envir = .GlobalEnv), add = TRUE)
  } else {
    on.exit(rm(".Random.seed", envir = .GlobalEnv, inherits = FALSE),
            add = TRUE)
  }
  set.seed(seed)

  if (!is.numeric(phi) || length(phi) != 1L || is.na(phi) || phi < 0) {
    stop("phi must be a single nonnegative numeric value.", call. = FALSE)
  }

  # ---- Default covariate spec for backward compatibility ----
  if (is.null(covariates)) {
    covariates <- list(
      Z = list(type = "binary", prob = 0.5),
      X = list(type = "continuous", dist = "normal", mean = 0, sd = 1)
    )
  }
  if (!is.list(covariates) || is.null(names(covariates)) ||
      any(!nzchar(names(covariates)))) {
    stop("`covariates` must be a named list, one entry per covariate.",
         call. = FALSE)
  }
  cov_names <- names(covariates)

  # ============================================================
  # OU helper (unchanged from previous version)
  # ============================================================
  generate_OU <- function(times, phi) {
    n_times <- length(times)
    if (n_times == 0) return(numeric(0))
    if (phi == 0) return(rep(0, n_times))
    if (n_times == 1) return(rnorm(1, 0, 1))
    ord <- order(times)
    times_sorted <- times[ord]
    xi <- numeric(n_times)
    xi[1] <- rnorm(1, 0, 1)
    for (j in 2:n_times) {
      dt  <- times_sorted[j] - times_sorted[j - 1]
      rho <- exp(-phi * dt)
      xi[j] <- rho * xi[j - 1] + sqrt(1 - rho^2) * rnorm(1)
    }
    xi[order(ord)]
  }

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
  # 1. Generate baseline covariates from the spec
  # ============================================================
  base_df <- data.frame(id = seq_len(n))
  for (cn in cov_names) {
    spec <- covariates[[cn]]
    if (is.null(spec$type) ||
        !spec$type %in% c("continuous", "binary")) {
      stop(sprintf("Covariate '%s' must declare type 'continuous' or 'binary'.",
                   cn), call. = FALSE)
    }
    if (spec$type == "binary") {
      p <- if (is.null(spec$prob)) 0.5 else spec$prob
      base_df[[cn]] <- rbinom(n, 1, p)
    } else {
      dist <- if (is.null(spec$dist)) "normal" else spec$dist
      if (dist == "normal") {
        m  <- if (is.null(spec$mean)) 0 else spec$mean
        sd <- if (is.null(spec$sd))   1 else spec$sd
        x <- rnorm(n, m, sd)
      } else if (dist == "uniform") {
        lo <- if (is.null(spec$min)) 0 else spec$min
        hi <- if (is.null(spec$max)) 1 else spec$max
        x <- runif(n, lo, hi)
      } else {
        stop(sprintf("Unknown dist '%s' for covariate '%s'.", dist, cn),
             call. = FALSE)
      }
      if (isTRUE(spec$standardize)) {
        # Keep the raw (pre-standardization) values alongside the standardized
        # column so descriptive summaries can recover original units. Stored
        # as `<name>_raw` and excluded from the design matrices used by the
        # data-generating model — those use only the standardized column.
        base_df[[paste0(cn, "_raw")]] <- x
        x <- as.numeric(scale(x))
      }
      base_df[[cn]] <- x
    }
  }

  # ============================================================
  # 2. Latent factor U and frailty eta
  # ============================================================
  U_i <- rnorm(n, 0, 1)
  mu0 <- -0.5 * sigma_zeta^2
  eta_i <- exp(mu0 + sigma_zeta * U_i)

  # ============================================================
  # 3. Linear-predictor helper: takes a named coefficient vector and
  #    returns the fixed-effect contribution sum_j coef_j * cov_j (plus
  #    intercept if present). Coefficient names not matching any
  #    declared covariate (other than 'intercept') trigger an error so
  #    typos surface immediately.
  # ============================================================
  fixed_lp <- function(coef, where) {
    if (is.null(names(coef)) || any(!nzchar(names(coef)))) {
      stop(sprintf("All entries of `%s` must be named.", where), call. = FALSE)
    }
    nm <- names(coef)
    bad <- setdiff(setdiff(nm, "intercept"), cov_names)
    if (length(bad) > 0L) {
      stop(sprintf("`%s` references covariates not in the spec: %s",
                   where, paste(bad, collapse = ", ")), call. = FALSE)
    }
    lp <- numeric(n)
    if ("intercept" %in% nm) lp <- lp + coef[["intercept"]]
    for (cn in setdiff(nm, "intercept")) {
      lp <- lp + coef[[cn]] * base_df[[cn]]
    }
    lp
  }

  # ============================================================
  # 4. Visit times
  # ============================================================
  tau <- time.end - time.start
  lambda_i <- eta_i * exp(fixed_lp(gamma, "gamma"))
  visit_times_list <- lapply(lambda_i, function(lam) {
    r_poisson_times(lam, tau) + time.start
  })
  n_visits <- vapply(visit_times_list, length, integer(1))

  # ============================================================
  # 5. Observation-process random effects q_i
  #    Random-effect design is [intercept, obs_random_covs].
  # ============================================================
  q_design_names <- c("intercept", obs_random_covs)
  for (req in c("delta", "sigma_q")) {
    v <- get(req)
    miss <- setdiff(q_design_names, names(v))
    if (length(miss) > 0L) {
      stop(sprintf("`%s` is missing entries for: %s", req,
                   paste(miss, collapse = ", ")), call. = FALSE)
    }
  }
  delta_vec   <- delta[q_design_names]
  sigma_q_vec <- sigma_q[q_design_names]
  q_dim <- length(q_design_names)
  q_resid <- matrix(rnorm(n * q_dim), nrow = n, ncol = q_dim)
  q_resid <- sweep(q_resid, 2, sigma_q_vec, "*")
  q_i <- outer(U_i, delta_vec) + q_resid
  colnames(q_i) <- paste0("q_", q_design_names)

  # ============================================================
  # 6. Outcome random effects b_i
  # ============================================================
  b_design_names <- c("intercept", outcome_random_covs)
  for (req in c("theta")) {
    v <- get(req)
    miss <- setdiff(b_design_names, names(v))
    if (length(miss) > 0L) {
      stop(sprintf("`%s` is missing entries for: %s", req,
                   paste(miss, collapse = ", ")), call. = FALSE)
    }
  }
  if (!is.matrix(Sigma_b) ||
      nrow(Sigma_b) != length(b_design_names) ||
      ncol(Sigma_b) != length(b_design_names)) {
    stop(sprintf("`Sigma_b` must be a %d x %d matrix.",
                 length(b_design_names), length(b_design_names)),
         call. = FALSE)
  }
  theta_vec <- theta[b_design_names]
  if (!requireNamespace("MASS", quietly = TRUE)) {
    stop("Package 'MASS' is required for sim_data_gen() (multivariate normal).",
         call. = FALSE)
  }
  b_resid <- MASS::mvrnorm(n, mu = rep(0, length(b_design_names)),
                           Sigma = Sigma_b)
  if (length(b_design_names) == 1L) b_resid <- matrix(b_resid, ncol = 1L)
  b_i <- outer(U_i, theta_vec) + b_resid
  colnames(b_i) <- paste0("b_", b_design_names)

  # ============================================================
  # 7. Pre-compute the fixed-effect contributions to the observation
  #    and outcome linear predictors at the subject level (they only
  #    depend on baseline covariates, so they don't change across visits
  #    for a given subject).
  # ============================================================
  alpha_lp <- fixed_lp(alpha, "alpha")  # length n
  beta_lp  <- fixed_lp(beta,  "beta")   # length n

  # Subject-level random-effect design matrices for the per-visit loop:
  # row i contains [1, base_df[i, obs_random_covs]]
  q_design_mat <- cbind(1, as.matrix(base_df[, obs_random_covs, drop = FALSE]))
  colnames(q_design_mat) <- q_design_names
  b_design_mat <- cbind(1, as.matrix(base_df[, outcome_random_covs, drop = FALSE]))
  colnames(b_design_mat) <- b_design_names

  # ============================================================
  # 8. Expand to long format
  # ============================================================
  # cov_names_long includes both the modeled covariate columns and any
  # `<name>_raw` columns generated alongside them (when `standardize = TRUE`).
  cov_names_long <- intersect(c(cov_names, paste0(cov_names, "_raw")),
                              names(base_df))
  long_list <- vector("list", n)
  for (i in seq_len(n)) {
    times_i <- visit_times_list[[i]]
    m_i <- length(times_i)

    # Per-subject covariate row, repeated m_i (or 1) times.
    # Important: reset row names before putting repeated covariate rows into
    # data.frame(); otherwise base R may warn about duplicated row names.
    cov_row <- base_df[i, cov_names_long, drop = FALSE]
    rownames(cov_row) <- NULL

    if (m_i == 0L) {
      long_list[[i]] <- data.frame(
        id   = i,
        cov_row,
        time = NA_real_,
        R    = NA_integer_,
        Y    = NA_real_,
        check.names = FALSE
      )
      rownames(long_list[[i]]) <- NULL
    } else {
      xi_i <- generate_OU(times_i, phi)

      # Observation: alpha_lp[i] + q_i %*% z_obs(i) + xi_i
      lin_random_obs <- sum(q_i[i, ] * q_design_mat[i, ])
      R_tilde <- alpha_lp[i] + lin_random_obs + xi_i
      R_i <- as.integer(R_tilde > 0)

      # Outcome: beta_lp[i] + beta_t * t + b_i %*% z_out(i) + epsilon
      lin_random_out <- sum(b_i[i, ] * b_design_mat[i, ])
      mu_Y <- beta_lp[i] + beta_t * times_i + lin_random_out
      Y_i <- mu_Y + rnorm(m_i, 0, sigma_y)
      Y_i[R_i == 0L] <- NA_real_

      cov_rep <- cov_row[rep(1L, m_i), , drop = FALSE]
      rownames(cov_rep) <- NULL

      long_list[[i]] <- data.frame(
        id   = rep(i, m_i),
        cov_rep,
        time = times_i,
        R    = R_i,
        Y    = Y_i,
        check.names = FALSE
      )
      rownames(long_list[[i]]) <- NULL
    }
  }

  # data.table::rbindlist is faster and ignores row names; fall back to base R
  # if data.table is not installed.
  if (requireNamespace("data.table", quietly = TRUE)) {
    long_df <- as.data.frame(data.table::rbindlist(long_list, use.names = TRUE, fill = TRUE))
  } else {
    long_df <- do.call(rbind, long_list)
  }
  rownames(long_df) <- NULL

  # Constant censoring time C, one column at end of frame
  long_df$C <- time.end

  # ============================================================
  # 9. Latent diagnostics frame (subject-level, not merged into long_df)
  # ============================================================
  latent_df <- data.frame(
    id = base_df$id,
    U  = U_i,
    eta = eta_i,
    n_visits = n_visits
  )
  latent_df <- cbind(latent_df, q_i, b_i)

  # ============================================================
  # 10. Summaries
  # ============================================================
  prop_no_visit <- mean(n_visits == 0)
  obs_rate <- mean(long_df$R[!is.na(long_df$R)])

  if (isTRUE(verbose)) {
    message(sprintf("Simulation summary (n=%d):", n))
    message(sprintf("  Subjects with no visits: %.1f%%", 100 * prop_no_visit))
    message(sprintf("  Total visits: %d", sum(n_visits)))
    message(sprintf("  Median visits per subject: %.1f", median(n_visits)))
    message(sprintf("  Observation rate: %.1f%%", 100 * obs_rate))
    message(sprintf("  Total observed outcomes: %d",
                    sum(long_df$R == 1L, na.rm = TRUE)))
  }

  structure(
    list(
      long_data     = long_df,
      latent        = latent_df,
      prop_no_visit = prop_no_visit,
      obs_rate      = obs_rate,
      params = list(
        n          = n,
        time.start = time.start,
        time.end   = time.end,
        seed       = seed,
        covariates = covariates,
        gamma      = gamma,
        sigma_zeta = sigma_zeta,
        alpha      = alpha,
        obs_random_covs = obs_random_covs,
        delta      = delta,
        sigma_q    = sigma_q,
        phi        = phi,
        beta       = beta,
        beta_t     = beta_t,
        outcome_random_covs = outcome_random_covs,
        theta      = theta,
        Sigma_b    = Sigma_b,
        sigma_y    = sigma_y
      )
    ),
    class = "sim_CIMEHR"
  )
}

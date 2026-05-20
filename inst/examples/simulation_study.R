## ---- simulations/simulation_study.R ------------------------------------------
## Monte Carlo simulation study evaluating bias, empirical SE, mean SE,
## SE ratio, RMSE, and 95% CI coverage for four methods on simulated
## Electronic Health Record (EHR) data. The design crosses two Z types
## (continuous = NSES, binary = Gender) with four values of the
## Ornstein-Uhlenbeck (OU) serial-correlation parameter (phi = 0, 0.001,
## 0.1, 0.3), producing eight scenarios. Each scenario runs N_SIM Monte
## Carlo replications; each replication generates a dataset of N_SUBJ
## participants and fits all selected methods.
##
## Data-generating parameters follow the collaborator-specified setting:
##   - The focal covariate (Z) takes the collaborator's Z values
##     (gamma_Z=0.3, alpha_Z=-0.6, beta_Z=-0.5, delta_Z=-0.3,
##     sigma_q_Z=0.2, theta_Z=0.0).
##   - Age plays the role of X with the collaborator's X values
##     (gamma_X=0.5, alpha_X=-0.3, beta_X=0.5).
##   - The four other covariates (whichever of NSES/Gender is not focal,
##     plus Marital, Black, Hispanic) retain their canonical "background"
##     values from the existing cohort calibration.
##   - Scalars: sigma_zeta=0.5, sigma_y=0.7, time.end=120, beta_t=0.01,
##     Sigma_b=diag(0.5, 0.8). Per-replicate seed: base_seed + 1e6 + sim
##     with base_seed = 20260508.
##
## Methods compared:
##   - CIMEHR                                            (this package)
##   - Joint_modeling_visiting_and_longitudinal_Liang    (this package, as CIMEHR_Liang)
##   - Inverse_intensity_rate_ratio_weighting            (this package, as CIMEHR_IIRR)
##   - IrregLong::iiwgee                                 (external comparator, as IrregLong_iiwgee)
##
## SE / coverage availability:
##   - CIMEHR             : analytic SE -> coverage available
##   - CIMEHR_Liang       : analytic SE -> coverage available
##   - CIMEHR_IIRR        : no analytic SE -> coverage reported as ---
##   - IrregLong_iiwgee   : analytic SE from GEE summary -> coverage available
##
## Tracked parameters:
##   beta on each of: Age, Gender, Marital, Black, Hispanic, NSES.
##   The intercept is excluded. True beta is focal-dependent: the focal
##   covariate carries the collaborator's beta_Z=-0.5; Age carries
##   beta_X=0.5; the other four carry their background betas.
##
## Usage:
##   Rscript simulations/simulation_study.R
##   Rscript simulations/simulation_study.R --cores 8 --nsim 100
##   Rscript simulations/simulation_study.R --scenarios continuous_phi0,binary_phi0.3
##   Rscript simulations/simulation_study.R --methods CIMEHR,CIMEHR_Liang,CIMEHR_IIRR,IrregLong_iiwgee
## ---------------------------------------------------------------------------

# ============================================================================
# 0. Setup
# ============================================================================

suppressPackageStartupMessages({
  library(CIMEHR)
  library(IrregLong)
  library(survival)
  library(geepack)
})

args      <- commandArgs(trailingOnly = TRUE)
N_SIM     <- 1000L
N_SUBJ    <- 2000L
N_CORES   <- 1L
OUT_FILE  <- "simulations/sim_results.rds"
TIME_END  <- 120
BASE_SEED <- 20260508L

ALL_SCENARIOS <- c(
  "continuous_phi0",     "continuous_phi0.001",
  "continuous_phi0.1",   "continuous_phi0.3",
  "binary_phi0",         "binary_phi0.001",
  "binary_phi0.1",       "binary_phi0.3"
)

ALL_METHODS <- c("CIMEHR", "CIMEHR_Liang", "CIMEHR_IIRR", "IrregLong_iiwgee")

SELECTED_SCENARIOS <- ALL_SCENARIOS
SELECTED_METHODS   <- ALL_METHODS

parse_arg <- function(flag) {
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) NULL else args[i + 1L]
}

if (!is.null(parse_arg("--cores"))) N_CORES  <- as.integer(parse_arg("--cores"))
if (!is.null(parse_arg("--nsim")))  N_SIM    <- as.integer(parse_arg("--nsim"))
if (!is.null(parse_arg("--n")))     N_SUBJ   <- as.integer(parse_arg("--n"))
if (!is.null(parse_arg("--out")))   OUT_FILE <- parse_arg("--out")

if (!is.null(parse_arg("--scenarios"))) {
  SELECTED_SCENARIOS <- strsplit(parse_arg("--scenarios"), ",", fixed = TRUE)[[1]]
  bad <- setdiff(SELECTED_SCENARIOS, ALL_SCENARIOS)
  if (length(bad) > 0L) stop("Unknown scenario(s): ", paste(bad, collapse = ", "))
}

if (!is.null(parse_arg("--methods"))) {
  SELECTED_METHODS <- strsplit(parse_arg("--methods"), ",", fixed = TRUE)[[1]]
  bad <- setdiff(SELECTED_METHODS, ALL_METHODS)
  if (length(bad) > 0L) stop("Unknown method(s): ", paste(bad, collapse = ", "))
}

cat(sprintf("Simulation study: %d replications, n = %d, cores = %d\n",
            N_SIM, N_SUBJ, N_CORES))
cat(sprintf("Scenarios: %s\n", paste(SELECTED_SCENARIOS, collapse = ", ")))
cat(sprintf("Methods:   %s\n", paste(SELECTED_METHODS, collapse = ", ")))
cat(sprintf("Output:    %s\n\n", OUT_FILE))


# ============================================================================
# 1. Data-generating parameters
# ============================================================================

EHR_COVS <- list(
  Age      = list(type = "continuous", dist = "uniform",
                  min = 18, max = 99, standardize = TRUE),
  Gender   = list(type = "binary", prob = 0.5),
  Marital  = list(type = "binary", prob = 0.4),
  Black    = list(type = "binary", prob = 0.30),
  Hispanic = list(type = "binary", prob = 0.20),
  NSES     = list(type = "continuous", dist = "normal",
                  mean = 0, sd = 1, standardize = FALSE)
)

GAMMA_Z   <- 0.3;   GAMMA_X   <- 0.5
ALPHA_Z   <- -0.6;  ALPHA_X   <- -0.3
BETA_Z    <- -0.5;  BETA_X    <- 0.5
DELTA_Z   <- -0.3
SIGMA_Q_Z <- 0.2
THETA_Z   <- 0.0

GAMMA_INTERCEPT   <- -2.8
ALPHA_INTERCEPT   <-  0.2
BETA_INTERCEPT    <- -2.0
DELTA_INTERCEPT   <- -0.3
SIGMA_Q_INTERCEPT <-  0.5
THETA_INTERCEPT   <-  0.4
SIGMA_ZETA        <-  0.5
SIGMA_Y           <-  0.7
BETA_T            <-  0.01
SIGMA_B           <- matrix(c(0.5, 0.0, 0.0, 0.8), 2, 2)

BACKGROUND_GAMMA <- c(
  Gender   = -0.10,
  Marital  = -0.05,
  Black    =  0.15,
  Hispanic = -0.05,
  NSES     =  0.12
)

BACKGROUND_ALPHA <- c(
  Gender   = -0.10,
  Marital  =  0.00,
  Black    =  0.05,
  Hispanic = -0.05,
  NSES     =  0.20
)

BACKGROUND_BETA <- c(
  Gender   =  0.10,
  Marital  = -0.03,
  Black    =  0.15,
  Hispanic =  0.05,
  NSES     = -0.25
)

PARAM_NAMES <- c("Age", "Gender", "Marital", "Black", "Hispanic", "NSES")
ALL_COVARS  <- PARAM_NAMES

build_true_params <- function(focal) {
  background_others <- setdiff(names(BACKGROUND_GAMMA), focal)

  gamma <- c(
    intercept = GAMMA_INTERCEPT,
    Age = GAMMA_X,
    stats::setNames(GAMMA_Z, focal),
    BACKGROUND_GAMMA[background_others]
  )

  alpha <- c(
    intercept = ALPHA_INTERCEPT,
    Age = ALPHA_X,
    stats::setNames(ALPHA_Z, focal),
    BACKGROUND_ALPHA[background_others]
  )

  beta <- c(
    intercept = BETA_INTERCEPT,
    Age = BETA_X,
    stats::setNames(BETA_Z, focal),
    BACKGROUND_BETA[background_others]
  )

  gamma <- c(gamma["intercept"], gamma[PARAM_NAMES])
  alpha <- c(alpha["intercept"], alpha[PARAM_NAMES])
  beta  <- c(beta["intercept"],  beta[PARAM_NAMES])

  list(gamma = gamma, alpha = alpha, beta = beta)
}


# ============================================================================
# 2. Method definitions
# ============================================================================

z_crit <- stats::qnorm(0.975)

prep_irreglong_iiwgee <- function(dat) {
  obs <- dat[!is.na(dat$R) & dat$R == 1L, , drop = FALSE]
  obs <- obs[order(obs$id, obs$time), , drop = FALSE]

  # Keep subjects with at least two observed visits.
  n_obs <- ave(obs$time, obs$id, FUN = length)
  obs <- obs[n_obs >= 2L, , drop = FALSE]

  # Remap to dense integer IDs while preserving repeated rows per subject.
  obs$id <- as.integer(factor(obs$id, levels = unique(obs$id)))
  rownames(obs) <- NULL

  obs$event <- 1L

  maxfu <- data.frame(
    id    = sort(unique(obs$id)),
    maxfu = TIME_END
  )

  list(obs = obs, maxfu = maxfu)
}

extract_geeglm_coef <- function(fit, params = PARAM_NAMES) {
  s <- summary(fit)
  coef_tab <- s$coefficients

  keep_rows <- intersect(params, rownames(coef_tab))

  estimate_col <- if ("Estimate" %in% colnames(coef_tab)) {
    "Estimate"
  } else {
    colnames(coef_tab)[1L]
  }

  se_col <- intersect(
    c("Std.err", "Std. Error", "Robust S.E.", "Robust SE", "Naive S.E."),
    colnames(coef_tab)
  )

  if (length(se_col) == 0L) {
    stop("Could not identify standard-error column in IrregLong::iiwgee GEE summary.")
  }

  se_col <- se_col[1L]

  beta <- stats::setNames(coef_tab[keep_rows, estimate_col], keep_rows)
  se   <- stats::setNames(coef_tab[keep_rows, se_col], keep_rows)

  list(beta = beta, se = se)
}

methods <- list(

  CIMEHR = function(dat, focal) {
    fit <- CIMEHR(
      data                          = dat,
      y_col                         = "log_HbA1c",
      covars_visit_XV               = ALL_COVARS,
      covars_outcome_fixed_XY       = ALL_COVARS,
      covars_outcome_random_link_ZY = focal,
      covars_obs_fixed_XO           = ALL_COVARS,
      covars_obs_random_link_ZO     = focal
    )

    beta <- fit$longitudinal_outcome$beta_Y
    se   <- fit$longitudinal_outcome$se_beta_naive

    keep <- intersect(PARAM_NAMES, names(beta))
    beta <- beta[keep]
    se   <- se[keep]

    list(
      beta    = beta,
      se      = se,
      ci_low  = beta - z_crit * se,
      ci_high = beta + z_crit * se,
      extra   = list(
        gamma = fit$visiting_process$gamma_V,
        sigma2_zeta = fit$visiting_process$sigma2_zeta
      )
    )
  },

  CIMEHR_Liang = function(dat, focal) {
    fit <- Joint_modeling_visiting_and_longitudinal_Liang(
      data        = dat,
      outcome     = "log_HbA1c",
      visit_covs  = ALL_COVARS,
      long_covs   = ALL_COVARS,
      random_covs = focal
    )

    beta <- fit$beta_hat
    keep <- intersect(PARAM_NAMES, names(beta))
    beta <- beta[keep]
    se   <- fit$se_beta[keep]

    list(
      beta    = beta,
      se      = se,
      ci_low  = beta - z_crit * se,
      ci_high = beta + z_crit * se
    )
  },

  CIMEHR_IIRR = function(dat, focal) {
    out_form <- stats::reformulate(
      c(ALL_COVARS, "time"),
      response = "log_HbA1c"
    )

    fit <- Inverse_intensity_rate_ratio_weighting(
      data       = dat,
      visit_covs = ALL_COVARS,
      formula    = out_form
    )

    beta <- fit$beta_hat
    keep <- intersect(PARAM_NAMES, names(beta))

    list(
      beta    = beta[keep],
      se      = NULL,
      ci_low  = NULL,
      ci_high = NULL
    )
  },

  IrregLong_iiwgee = function(dat, focal) {
    prep <- prep_irreglong_iiwgee(dat)

    out_form <- stats::reformulate(
      c(ALL_COVARS, "time"),
      response = "log_HbA1c"
    )

    visit_form <- survival::Surv(time.lag, time, event) ~
      Age + Gender + Marital + Black + Hispanic + NSES + cluster(id)

    fit <- suppressWarnings(
      IrregLong::iiwgee(
        formulagee = out_form,
        formulaph  = visit_form,
        data       = prep$obs,
        id         = "id",
        time       = "time",
        event      = "event",
        invariant  = c("id", ALL_COVARS),
        lagvars    = "time",
        lagfirst   = 0,
        maxfu      = prep$maxfu,
        first      = TRUE
      )
    )

    coef_out <- extract_geeglm_coef(fit$geefit, PARAM_NAMES)

    beta <- coef_out$beta
    se   <- coef_out$se

    list(
      beta    = beta,
      se      = se,
      ci_low  = beta - z_crit * se,
      ci_high = beta + z_crit * se
    )
  }
)


# ============================================================================
# 3. Single-replication worker
# ============================================================================

make_scenario_params <- function(focal) {
  if (!focal %in% ALL_COVARS) {
    stop("`focal` must be one of: ", paste(ALL_COVARS, collapse = ", "),
         call. = FALSE)
  }

  tp <- build_true_params(focal)

  list(
    gamma   = tp$gamma,
    alpha   = tp$alpha,
    beta    = tp$beta,
    obs_random_covs     = focal,
    outcome_random_covs = focal,
    delta   = c(intercept = DELTA_INTERCEPT,
                stats::setNames(DELTA_Z, focal)),
    sigma_q = c(intercept = SIGMA_Q_INTERCEPT,
                stats::setNames(SIGMA_Q_Z, focal)),
    theta   = c(intercept = THETA_INTERCEPT,
                stats::setNames(THETA_Z, focal)),
    Sigma_b = SIGMA_B
  )
}

run_one_rep <- function(sim_id, Z_type, phi, focal) {
  sp <- make_scenario_params(focal)

  rep_seed <- BASE_SEED + 1000000L + sim_id

  sim <- suppressMessages(
    sim_data_gen(
      n          = N_SUBJ,
      time.start = 0,
      time.end   = TIME_END,
      seed       = rep_seed,
      covariates = EHR_COVS,
      gamma      = sp$gamma,
      sigma_zeta = SIGMA_ZETA,
      alpha      = sp$alpha,
      obs_random_covs = sp$obs_random_covs,
      delta      = sp$delta,
      sigma_q    = sp$sigma_q,
      phi        = phi,
      beta       = sp$beta,
      beta_t     = BETA_T,
      outcome_random_covs = sp$outcome_random_covs,
      theta      = sp$theta,
      Sigma_b    = sp$Sigma_b,
      sigma_y    = SIGMA_Y,
      verbose    = FALSE
    )
  )

  dat <- sim$long_data
  rownames(dat) <- NULL

  dat$HbA1c     <- exp(log(5.6) + 0.05 * dat$Y)
  dat$log_HbA1c <- log(dat$HbA1c)

  rep_results <- list()

  for (m in SELECTED_METHODS) {
    rep_results[[m]] <- tryCatch(
      methods[[m]](dat, focal = focal),
      error = function(e) {
        list(
          beta    = NULL,
          se      = NULL,
          ci_low  = NULL,
          ci_high = NULL,
          error   = conditionMessage(e)
        )
      }
    )
  }

  list(sim_id = sim_id, results = rep_results)
}


# ============================================================================
# 4. Run selected scenarios
# ============================================================================

scenario_spec <- function(s) {
  switch(
    s,
    "continuous_phi0"     = list(Z_type = "continuous", phi = 0,     focal = "NSES"),
    "continuous_phi0.001" = list(Z_type = "continuous", phi = 0.001, focal = "NSES"),
    "continuous_phi0.1"   = list(Z_type = "continuous", phi = 0.1,   focal = "NSES"),
    "continuous_phi0.3"   = list(Z_type = "continuous", phi = 0.3,   focal = "NSES"),
    "binary_phi0"         = list(Z_type = "binary",     phi = 0,     focal = "Gender"),
    "binary_phi0.001"     = list(Z_type = "binary",     phi = 0.001, focal = "Gender"),
    "binary_phi0.1"       = list(Z_type = "binary",     phi = 0.1,   focal = "Gender"),
    "binary_phi0.3"       = list(Z_type = "binary",     phi = 0.3,   focal = "Gender"),
    stop("Unknown scenario: ", s)
  )
}

run_scenario <- function(scenario_label) {
  spec <- scenario_spec(scenario_label)

  cat(sprintf(
    "\nRunning scenario: %s (Z_type = %s, phi = %.3f) ...\n",
    scenario_label, spec$Z_type, spec$phi
  ))

  if (N_CORES > 1L && requireNamespace("parallel", quietly = TRUE)) {
    cat(sprintf("  Using %d cores\n", N_CORES))

    reps <- parallel::mclapply(
      seq_len(N_SIM),
      function(i) {
        run_one_rep(
          sim_id = i,
          Z_type = spec$Z_type,
          phi    = spec$phi,
          focal  = spec$focal
        )
      },
      mc.cores = N_CORES,
      mc.set.seed = TRUE
    )

  } else {
    reps <- vector("list", N_SIM)

    for (i in seq_len(N_SIM)) {
      cat(sprintf("  Replicate %d / %d\r", i, N_SIM))
      flush.console()

      reps[[i]] <- run_one_rep(
        sim_id = i,
        Z_type = spec$Z_type,
        phi    = spec$phi,
        focal  = spec$focal
      )
    }

    cat(sprintf("  Replicate %d / %d -- done\n", N_SIM, N_SIM))
  }

  list(scenario = scenario_label, spec = spec, reps = reps)
}

cat("Starting simulations ...\n")
t0 <- Sys.time()

all_scenario_runs <- lapply(SELECTED_SCENARIOS, run_scenario)

elapsed <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
cat(sprintf("\nAll scenarios done in %.1f minutes.\n\n", elapsed))


# ============================================================================
# 5. Assemble results
# ============================================================================

assemble <- function(scenario_run) {
  scen   <- scenario_run$scenario
  Z_type <- scenario_run$spec$Z_type
  phi    <- scenario_run$spec$phi
  focal  <- scenario_run$spec$focal

  rows <- list()

  for (rep in scenario_run$reps) {
    for (m in SELECTED_METHODS) {
      res <- rep$results[[m]]
      failed <- is.null(res$beta)

      for (p in PARAM_NAMES) {
        rows[[length(rows) + 1L]] <- data.frame(
          scenario = scen,
          Z_type   = Z_type,
          phi      = phi,
          focal    = focal,
          sim_id   = rep$sim_id,
          method   = m,
          param    = p,
          est      = if (!failed && p %in% names(res$beta)) {
            res$beta[[p]]
          } else {
            NA_real_
          },
          se       = if (!failed && !is.null(res$se) && p %in% names(res$se)) {
            res$se[[p]]
          } else {
            NA_real_
          },
          ci_low   = if (!failed && !is.null(res$ci_low) && p %in% names(res$ci_low)) {
            res$ci_low[[p]]
          } else {
            NA_real_
          },
          ci_high  = if (!failed && !is.null(res$ci_high) && p %in% names(res$ci_high)) {
            res$ci_high[[p]]
          } else {
            NA_real_
          },
          failed   = failed,
          stringsAsFactors = FALSE
        )
      }
    }
  }

  do.call(rbind, rows)
}

sim_df <- do.call(rbind, lapply(all_scenario_runs, assemble))
rownames(sim_df) <- NULL

.true_beta_by_scenario <- function(scen) {
  spec <- scenario_spec(scen)
  build_true_params(spec$focal)$beta
}

.true_beta_lookup <- lapply(unique(sim_df$scenario), function(s) {
  tb <- .true_beta_by_scenario(s)
  stats::setNames(as.numeric(tb[PARAM_NAMES]), PARAM_NAMES)
})

names(.true_beta_lookup) <- unique(sim_df$scenario)

sim_df$true_value <- vapply(
  seq_len(nrow(sim_df)),
  function(i) .true_beta_lookup[[sim_df$scenario[i]]][[sim_df$param[i]]],
  numeric(1)
)


# ============================================================================
# 6. Compute performance metrics
# ============================================================================

compute_metrics <- function(df) {
  n_fail <- sum(df$failed, na.rm = TRUE)
  df     <- df[!df$failed & !is.na(df$est), , drop = FALSE]
  n_ok   <- nrow(df)

  if (n_ok == 0L) {
    return(data.frame(
      n_reps   = 0L,
      n_fail   = n_fail,
      bias     = NA_real_,
      emp_se   = NA_real_,
      mean_se  = NA_real_,
      se_ratio = NA_real_,
      rmse     = NA_real_,
      coverage = NA_real_
    ))
  }

  true   <- df$true_value[1L]
  bias   <- mean(df$est) - true
  emp_se <- stats::sd(df$est)
  rmse   <- sqrt(mean((df$est - true)^2))

  has_se   <- !is.na(df$se)
  mean_se  <- if (any(has_se)) mean(df$se[has_se]) else NA_real_
  se_ratio <- if (!is.na(mean_se) && emp_se > 0) mean_se / emp_se else NA_real_

  has_ci <- !is.na(df$ci_low) & !is.na(df$ci_high)
  coverage <- if (any(has_ci)) {
    mean(df$ci_low[has_ci] <= true & true <= df$ci_high[has_ci])
  } else {
    NA_real_
  }

  data.frame(
    n_reps   = n_ok,
    n_fail   = n_fail,
    bias     = bias,
    emp_se   = emp_se,
    mean_se  = mean_se,
    se_ratio = se_ratio,
    rmse     = rmse,
    coverage = coverage
  )
}

metrics_keys <- paste(sim_df$scenario, sim_df$method, sim_df$param, sep = "||")

metrics <- do.call(
  rbind,
  lapply(split(sim_df, metrics_keys), compute_metrics)
)

keys <- strsplit(rownames(metrics), "\\|\\|")
metrics$scenario <- vapply(keys, `[`, character(1), 1)
metrics$method   <- vapply(keys, `[`, character(1), 2)
metrics$param    <- vapply(keys, `[`, character(1), 3)
rownames(metrics) <- NULL

spec_lookup <- do.call(
  rbind,
  lapply(SELECTED_SCENARIOS, function(s) {
    sp <- scenario_spec(s)
    data.frame(
      scenario = s,
      Z_type   = sp$Z_type,
      phi      = sp$phi,
      focal    = sp$focal,
      stringsAsFactors = FALSE
    )
  })
)

metrics <- merge(
  metrics,
  spec_lookup,
  by = "scenario",
  all.x = TRUE,
  sort = FALSE
)

metrics <- metrics[, c(
  "scenario", "Z_type", "phi", "focal",
  "method", "param", "n_reps", "n_fail",
  "bias", "emp_se", "mean_se",
  "se_ratio", "rmse", "coverage"
)]

metrics <- metrics[
  order(
    match(metrics$scenario, SELECTED_SCENARIOS),
    match(metrics$method, SELECTED_METHODS),
    match(metrics$param, PARAM_NAMES)
  ),
]


# ============================================================================
# 7. Print results
# ============================================================================

cat(strrep("=", 80), "\n", sep = "")
cat("SIMULATION RESULTS\n")
cat(sprintf("  Replications: %d | Subjects per rep: %d\n", N_SIM, N_SUBJ))
cat("  True beta varies by scenario: the focal covariate gets the collaborator's\n")
cat("  Z value; Age always gets the X value; the four remaining covariates\n")
cat("  keep their background values. See per-scenario blocks below.\n")
cat(strrep("=", 80), "\n\n", sep = "")

for (scen in SELECTED_SCENARIOS) {
  spec <- scenario_spec(scen)

  cat(sprintf(
    "--- %s (Z_type = %s, phi = %.3f, focal = %s) ---\n\n",
    scen, spec$Z_type, spec$phi, spec$focal
  ))

  scen_rows  <- metrics[metrics$scenario == scen, , drop = FALSE]
  true_betas <- .true_beta_lookup[[scen]]

  for (p in PARAM_NAMES) {
    sub <- scen_rows[scen_rows$param == p, , drop = FALSE]
    if (nrow(sub) == 0L) next

    is_focal <- p == spec$focal

    cat(sprintf(
      "  beta_%s (true = %+.3f)%s\n\n",
      p,
      true_betas[[p]],
      if (is_focal) "  [FOCAL]" else ""
    ))

    print_df <- data.frame(
      Method   = sub$method,
      Reps     = sub$n_reps,
      Fail     = sub$n_fail,
      Bias     = sprintf("%+.4f", sub$bias),
      Emp.SE   = sprintf("%.4f",  sub$emp_se),
      Mean.SE  = ifelse(is.na(sub$mean_se),  "---", sprintf("%.4f", sub$mean_se)),
      SE.Ratio = ifelse(is.na(sub$se_ratio), "---", sprintf("%.3f", sub$se_ratio)),
      RMSE     = sprintf("%.4f", sub$rmse),
      Cov.95   = ifelse(is.na(sub$coverage), "---", sprintf("%.3f", sub$coverage))
    )

    print(print_df, row.names = FALSE, right = FALSE)
    cat("\n")
  }
}

cat("Notes:\n")
cat("  - SE.Ratio = Mean.SE / Emp.SE; values near 1.0 suggest calibrated SEs.\n")
cat("  - Cov.95   = empirical coverage of the 95% CI; values near 0.950 are expected.\n")
cat("  - Mean.SE / SE.Ratio / Cov.95 are reported as --- for CIMEHR_IIRR,\n")
cat("    because analytic standard errors are not returned by that method.\n")
cat("    Bootstrap inference would be needed for uncertainty quantification.\n\n")


# ============================================================================
# 8. Save everything
# ============================================================================

dir.create(dirname(OUT_FILE), recursive = TRUE, showWarnings = FALSE)

true_params_by_scenario <- lapply(SELECTED_SCENARIOS, function(s) {
  spec <- scenario_spec(s)
  tp   <- build_true_params(spec$focal)

  list(
    focal      = spec$focal,
    phi        = spec$phi,
    gamma      = tp$gamma,
    alpha      = tp$alpha,
    beta       = tp$beta,
    delta      = c(
      intercept = DELTA_INTERCEPT,
      stats::setNames(DELTA_Z, spec$focal)
    ),
    sigma_q    = c(
      intercept = SIGMA_Q_INTERCEPT,
      stats::setNames(SIGMA_Q_Z, spec$focal)
    ),
    theta      = c(
      intercept = THETA_INTERCEPT,
      stats::setNames(THETA_Z, spec$focal)
    ),
    Sigma_b    = SIGMA_B,
    sigma_zeta = SIGMA_ZETA,
    sigma_y    = SIGMA_Y,
    beta_t     = BETA_T,
    time_end   = TIME_END
  )
})

names(true_params_by_scenario) <- SELECTED_SCENARIOS

saveRDS(
  list(
    sim_df  = sim_df,
    metrics = metrics,
    settings = list(
      n_sim     = N_SIM,
      n_subj    = N_SUBJ,
      base_seed = BASE_SEED,
      time_end  = TIME_END,
      scenarios = SELECTED_SCENARIOS,
      methods   = SELECTED_METHODS,
      anchors = list(
        gamma_Z = GAMMA_Z,
        gamma_X = GAMMA_X,
        gamma_intercept = GAMMA_INTERCEPT,
        alpha_Z = ALPHA_Z,
        alpha_X = ALPHA_X,
        alpha_intercept = ALPHA_INTERCEPT,
        beta_Z = BETA_Z,
        beta_X = BETA_X,
        beta_intercept = BETA_INTERCEPT,
        delta_Z = DELTA_Z,
        sigma_q_Z = SIGMA_Q_Z,
        theta_Z = THETA_Z,
        sigma_zeta = SIGMA_ZETA,
        sigma_y = SIGMA_Y,
        Sigma_b = SIGMA_B,
        beta_t = BETA_T
      ),
      background_gamma = BACKGROUND_GAMMA,
      background_alpha = BACKGROUND_ALPHA,
      background_beta  = BACKGROUND_BETA,
      true_params_by_scenario = true_params_by_scenario
    ),
    elapsed_min = elapsed
  ),
  file = OUT_FILE
)

cat(sprintf("Results saved to %s\n", OUT_FILE))

## ---- simulations/simulation_study.R ------------------------------------------
## Monte Carlo simulation study evaluating bias, empirical SE, mean SE,
## SE ratio, RMSE, and 95% CI coverage for five methods on simulated EHR
## data. The design crosses two Z types (continuous = NSES, binary = Gender)
## with two values of the OU serial-correlation parameter (phi = 0, 0.3),
## producing four scenarios. Each scenario runs N_SIM Monte Carlo
## replications; each replication generates a fresh dataset of N_SUBJ
## participants and fits all methods.
##
## Methods compared:
##   - CIMEHR                                            (this package)
##   - Joint_modeling_visiting_and_longitudinal_Liang    (this package)
##   - IrregLong::Liang                                  (IrregLong)
##   - Inverse_intensity_rate_ratio_weighting            (this package)
##   - IrregLong::iiwgee                                 (IrregLong)
##
## SE / coverage availability:
##   - CIMEHR              : sandwich SE -> coverage available
##   - IrregLong::iiwgee   : sandwich SE -> coverage available
##   - All others          : no analytic SE -> coverage reported as ---
##
## Tracked parameters:
##   beta on each of: Age, Gender, Marital, Black, Hispanic, NSES.
##   The intercept is excluded.
##
## Usage:
##   Rscript simulations/simulation_study.R                   # all four scenarios
##   Rscript simulations/simulation_study.R --cores 8 --nsim 100
##   Rscript simulations/simulation_study.R --scenarios continuous_phi0,binary_phi0.3
##   Rscript simulations/simulation_study.R --methods CIMEHR,IrregLong_iiwgee
## ---------------------------------------------------------------------------

# ============================================================================
# 0. Setup
# ============================================================================

suppressPackageStartupMessages({
  library(CIMEHR)
  library(IrregLong)
  library(survival)   # Surv, frailty, cluster
  library(geepack)    # iiwgee dependency
  library(data.table) # IrregLong dependency
})

args      <- commandArgs(trailingOnly = TRUE)
N_SIM     <- 1000L
N_SUBJ    <- 2000L
N_CORES   <- 1L
OUT_FILE  <- "simulations/sim_results.rds"
TIME_END  <- 60

ALL_SCENARIOS <- c("continuous_phi0", "continuous_phi0.3",
                   "binary_phi0",     "binary_phi0.3")
ALL_METHODS   <- c("CIMEHR", "CIMEHR_Liang", "IrregLong_Liang",
                   "CIMEHR_IIRR", "IrregLong_iiwgee")

SELECTED_SCENARIOS <- ALL_SCENARIOS
SELECTED_METHODS   <- ALL_METHODS

parse_arg <- function(flag) {
  i <- match(flag, args)
  if (is.na(i)) NULL else args[i + 1L]
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
cat(sprintf("Methods:   %s\n", paste(SELECTED_METHODS,   collapse = ", ")))
cat(sprintf("Output:    %s\n\n", OUT_FILE))

# ============================================================================
# 1. Data-generating parameters (true values for bias / coverage)
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

# True coefficients used by sim_data_gen() in every replication
TRUE_GAMMA <- c(intercept = -2.0, Age = 0.10, Gender = -0.10, Marital = -0.05,
                Black = 0.15, Hispanic = -0.05, NSES = 0.12)
TRUE_ALPHA <- c(intercept = 0.5,  Age = 0.10, Gender = -0.10, Marital = 0.00,
                Black = 0.05, Hispanic = -0.05, NSES = 0.20)
TRUE_BETA  <- c(intercept = -2.0, Age = 0.05, Gender = 0.10, Marital = -0.03,
                Black = 0.15, Hispanic = 0.05, NSES = -0.25)
TRUE_DELTA   <- c(intercept = -0.3, NSES = -0.5)
TRUE_SIGMA_Q <- c(intercept =  0.3, NSES =  0.2)
TRUE_THETA   <- c(intercept =  0.5, NSES =  0.30)
TRUE_SIGMA_B <- matrix(c(1, 0, 0, 2), 2, 2)

PARAM_NAMES <- c("Age", "Gender", "Marital", "Black", "Hispanic", "NSES")
ALL_COVARS  <- PARAM_NAMES                  # all six are model covariates


# ============================================================================
# 2. Method definitions
# ============================================================================
# Each method function takes a data.frame and returns a named list:
#   $beta    - named numeric vector of outcome coefficient estimates
#   $se      - named numeric vector of SE (NULL if unavailable)
#   $ci_low  - named numeric vector of CI lower bound (NULL if unavailable)
#   $ci_high - named numeric vector of CI upper bound (NULL if unavailable)
#   $extra   - optional list of additional estimates (e.g., gamma)

z_crit <- stats::qnorm(0.975)

# Helper: Surv-based formula for IrregLong
make_irreg_formula <- function(rhs_terms, lhs = "Surv(time.lag, time, event)",
                               include_frailty = FALSE,
                               include_cluster = FALSE) {
  rhs <- paste(rhs_terms, collapse = " + ")
  if (include_frailty) rhs <- paste0(rhs, " + frailty(id)")
  if (include_cluster) rhs <- paste0(rhs, " + cluster(id)")
  stats::as.formula(paste(lhs, "~", rhs))
}

# Build observed-data view for IrregLong (one row per visit with R == 1)
prep_irreglong <- function(data) {
  obs <- data[!is.na(data$R) & data$R == 1L, , drop = FALSE]
  obs <- obs[order(obs$id, obs$time), , drop = FALSE]
  obs$event <- 1L
  maxfu <- data.frame(id = unique(obs$id), maxfu = TIME_END)
  list(obs = obs, maxfu = maxfu)
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
      extra   = list(gamma = fit$visiting_process$gamma_V,
                     sigma2_zeta = fit$visiting_process$sigma2_zeta)
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
    list(beta = beta, se = NULL, ci_low = NULL, ci_high = NULL)
  },

  IrregLong_Liang = function(dat, focal) {
    prep <- prep_irreglong(dat)
    # IrregLong::Liang sets row.names from `id` internally and emits a
    # benign "non-unique row.names" warning when subjects have multiple
    # visits. The estimates are unaffected; suppress to keep stdout clean.
    fit <- suppressWarnings(IrregLong::Liang(
      data       = prep$obs,
      Yname      = "log_HbA1c",
      Xnames     = ALL_COVARS,
      Wnames     = focal,
      Znames     = ALL_COVARS,
      formulaobs = make_irreg_formula(ALL_COVARS, include_frailty = TRUE),
      id         = "id",
      time       = "time",
      invariant  = c("id", ALL_COVARS),
      lagvars    = "time",
      lagfirst   = 0,
      maxfu      = prep$maxfu,
      baseline   = 0
    ))
    b   <- as.numeric(fit)
    nm  <- names(fit)
    if (is.null(nm)) nm <- ALL_COVARS[seq_along(b)]
    names(b) <- nm
    keep <- intersect(PARAM_NAMES, names(b))
    list(beta = b[keep], se = NULL, ci_low = NULL, ci_high = NULL)
  },

  CIMEHR_IIRR = function(dat, focal) {
    out_form <- stats::reformulate(c(ALL_COVARS, "time"), response = "log_HbA1c")
    fit <- Inverse_intensity_rate_ratio_weighting(
      data       = dat,
      visit_covs = ALL_COVARS,
      formula    = out_form
    )
    beta <- fit$beta_hat
    keep <- intersect(PARAM_NAMES, names(beta))
    list(beta = beta[keep], se = NULL, ci_low = NULL, ci_high = NULL)
  },

  IrregLong_iiwgee = function(dat, focal) {
    prep        <- prep_irreglong(dat)
    out_form    <- stats::reformulate(c(ALL_COVARS, "time"), response = "log_HbA1c")
    visit_form  <- make_irreg_formula(ALL_COVARS, include_cluster = TRUE)
    # See note above on IrregLong row.names warning.
    fit <- suppressWarnings(IrregLong::iiwgee(
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
    ))
    s    <- summary(fit$geefit)
    coef_tab <- s$coefficients
    keep_rows <- intersect(PARAM_NAMES, rownames(coef_tab))
    beta <- setNames(coef_tab[keep_rows, "Estimate"], keep_rows)
    se   <- setNames(coef_tab[keep_rows, "Std.err"],  keep_rows)
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

# Generate one dataset and fit selected methods. Returns a list with:
#   $sim_id  - replicate id
#   $results - per-method list of $beta / $se / $ci_low / $ci_high / $error
make_scenario_params <- function(focal) {
  if (!focal %in% ALL_COVARS) {
    stop("`focal` must be one of: ", paste(ALL_COVARS, collapse = ", "),
         call. = FALSE)
  }
  list(
    obs_random_covs     = focal,
    outcome_random_covs = focal,
    delta   = c(intercept = -0.3, stats::setNames(-0.5,  focal)),
    sigma_q = c(intercept =  0.3, stats::setNames( 0.2,  focal)),
    theta   = c(intercept =  0.5, stats::setNames( 0.30, focal)),
    Sigma_b = matrix(c(1, 0, 0, 2), 2, 2)
  )
}

# Generate one dataset and fit selected methods. Returns a list with:
#   $sim_id  - replicate id
#   $results - per-method list of $beta / $se / $ci_low / $ci_high / $error
run_one_rep <- function(sim_id, Z_type, phi, focal) {
  sp <- make_scenario_params(focal)

  sim <- suppressMessages(
    sim_data_gen(
      n          = N_SUBJ,
      time.start = 0,
      time.end   = TIME_END,
      seed       = sim_id,
      covariates = EHR_COVS,
      gamma      = TRUE_GAMMA,
      sigma_zeta = 1.0,
      alpha      = TRUE_ALPHA,
      obs_random_covs = sp$obs_random_covs,
      delta      = sp$delta,
      sigma_q    = sp$sigma_q,
      phi        = phi,
      beta       = TRUE_BETA,
      beta_t     = 0.01,
      outcome_random_covs = sp$outcome_random_covs,
      theta      = sp$theta,
      Sigma_b    = sp$Sigma_b,
      sigma_y    = 1,
      verbose    = FALSE
    )
  )
  dat <- sim$long_data
  rownames(dat) <- NULL

  # Add raw and log-scale HbA1c columns to match the bundled `sim_ehr_data`.
  # `Y` is the standardized linear outcome; the back-transform mirrors
  # `data-raw/sim_ehr_data.R`. Methods in this script reference `log_HbA1c`.
  dat$HbA1c     <- exp(log(5.6) + 0.05 * dat$Y)
  dat$log_HbA1c <- log(dat$HbA1c)

  rep_results <- list()
  for (m in SELECTED_METHODS) {
    rep_results[[m]] <- tryCatch(
      methods[[m]](dat, focal = focal),
      error = function(e) {
        list(beta = NULL, se = NULL, ci_low = NULL, ci_high = NULL,
             error = conditionMessage(e))
      }
    )
  }
  list(sim_id = sim_id, results = rep_results)
}

# ============================================================================
# 4. Run all selected scenarios
# ============================================================================

# Scenario string -> (Z_type, phi, focal) lookup. The focal covariate is used
# in both the data-generating model and the fitted models: NSES for the
# continuous-Z scenarios and Gender for the binary-Z scenarios.
scenario_spec <- function(s) {
  switch(s,
         "continuous_phi0"   = list(Z_type = "continuous", phi = 0,   focal = "NSES"),
         "continuous_phi0.3" = list(Z_type = "continuous", phi = 0.3, focal = "NSES"),
         "binary_phi0"       = list(Z_type = "binary",     phi = 0,   focal = "Gender"),
         "binary_phi0.3"     = list(Z_type = "binary",     phi = 0.3, focal = "Gender"),
         stop("Unknown scenario: ", s)
  )
}

run_scenario <- function(scenario_label) {
  spec <- scenario_spec(scenario_label)
  cat(sprintf("\nRunning scenario: %s (Z_type = %s, phi = %.1f) ...\n",
              scenario_label, spec$Z_type, spec$phi))

  if (N_CORES > 1L && requireNamespace("parallel", quietly = TRUE)) {
    cat(sprintf("  Using %d cores\n", N_CORES))
    reps <- parallel::mclapply(
      seq_len(N_SIM),
      function(i) run_one_rep(i, Z_type = spec$Z_type, phi = spec$phi, focal = spec$focal),
      mc.cores = N_CORES, mc.set.seed = TRUE
    )
  } else {
    reps <- vector("list", N_SIM)
    for (i in seq_len(N_SIM)) {
      cat(sprintf("  Replicate %d / %d\r", i, N_SIM))
      flush.console()
      reps[[i]] <- run_one_rep(i, Z_type = spec$Z_type, phi = spec$phi, focal = spec$focal)
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
# 5. Assemble results into a tidy data.frame
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
          est      = if (!failed && p %in% names(res$beta))    res$beta[[p]]    else NA_real_,
          se       = if (!failed && !is.null(res$se)     && p %in% names(res$se))     res$se[[p]]     else NA_real_,
          ci_low   = if (!failed && !is.null(res$ci_low) && p %in% names(res$ci_low)) res$ci_low[[p]] else NA_real_,
          ci_high  = if (!failed && !is.null(res$ci_high)&& p %in% names(res$ci_high))res$ci_high[[p]]else NA_real_,
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
sim_df$true_value <- TRUE_BETA[sim_df$param]


# ============================================================================
# 6. Compute performance metrics per (scenario, method, param)
# ============================================================================

compute_metrics <- function(df) {
  n_fail <- sum(df$failed, na.rm = TRUE)
  df     <- df[!df$failed & !is.na(df$est), , drop = FALSE]
  n_ok   <- nrow(df)

  if (n_ok == 0L) {
    return(data.frame(
      n_reps = 0L, n_fail = n_fail,
      bias = NA_real_, emp_se = NA_real_, mean_se = NA_real_,
      se_ratio = NA_real_, rmse = NA_real_, coverage = NA_real_
    ))
  }

  true     <- df$true_value[1L]
  bias     <- mean(df$est) - true
  emp_se   <- stats::sd(df$est)
  rmse     <- sqrt(mean((df$est - true)^2))

  has_se   <- !is.na(df$se)
  mean_se  <- if (any(has_se)) mean(df$se[has_se]) else NA_real_
  se_ratio <- if (!is.na(mean_se) && emp_se > 0) mean_se / emp_se else NA_real_

  has_ci   <- !is.na(df$ci_low) & !is.na(df$ci_high)
  coverage <- if (any(has_ci)) {
    mean(df$ci_low[has_ci] <= true & true <= df$ci_high[has_ci])
  } else NA_real_

  data.frame(
    n_reps = n_ok, n_fail = n_fail,
    bias = bias, emp_se = emp_se, mean_se = mean_se,
    se_ratio = se_ratio, rmse = rmse, coverage = coverage
  )
}

metrics_keys <- paste(sim_df$scenario, sim_df$method, sim_df$param, sep = "||")
metrics <- do.call(rbind, lapply(split(sim_df, metrics_keys), compute_metrics))
keys <- strsplit(rownames(metrics), "\\|\\|")
metrics$scenario <- vapply(keys, `[`, character(1), 1)
metrics$method   <- vapply(keys, `[`, character(1), 2)
metrics$param    <- vapply(keys, `[`, character(1), 3)
rownames(metrics) <- NULL

# Attach Z_type, phi, focal columns for convenience
spec_lookup <- do.call(rbind, lapply(SELECTED_SCENARIOS, function(s) {
  sp <- scenario_spec(s)
  data.frame(scenario = s, Z_type = sp$Z_type, phi = sp$phi,
             focal = sp$focal, stringsAsFactors = FALSE)
}))
metrics <- merge(metrics, spec_lookup, by = "scenario", all.x = TRUE,
                 sort = FALSE)

metrics <- metrics[, c("scenario", "Z_type", "phi", "focal",
                       "method", "param", "n_reps", "n_fail",
                       "bias", "emp_se", "mean_se",
                       "se_ratio", "rmse", "coverage")]
metrics <- metrics[order(match(metrics$scenario, SELECTED_SCENARIOS),
                         match(metrics$method,   SELECTED_METHODS),
                         match(metrics$param,    PARAM_NAMES)), ]


# ============================================================================
# 7. Print results
# ============================================================================

cat(strrep("=", 80), "\n", sep = "")
cat("SIMULATION RESULTS\n")
cat(sprintf("  Replications: %d | Subjects per rep: %d\n", N_SIM, N_SUBJ))
cat("  True beta:\n")
for (p in PARAM_NAMES) {
  cat(sprintf("    %-9s = %+.3f\n", p, TRUE_BETA[[p]]))
}
cat(strrep("=", 80), "\n\n", sep = "")

for (scen in SELECTED_SCENARIOS) {
  spec <- scenario_spec(scen)
  cat(sprintf("--- %s (Z_type = %s, phi = %.1f, focal = %s) ---\n\n",
              scen, spec$Z_type, spec$phi, spec$focal))
  scen_rows <- metrics[metrics$scenario == scen, , drop = FALSE]

  for (p in PARAM_NAMES) {
    sub <- scen_rows[scen_rows$param == p, , drop = FALSE]
    if (nrow(sub) == 0L) next
    is_focal <- (p == spec$focal)
    cat(sprintf("  beta_%s (true = %+.3f)%s\n\n", p, TRUE_BETA[[p]],
                if (is_focal) "  [FOCAL]" else ""))

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
cat("  - SE.Ratio = Mean.SE / Emp.SE (should be near 1.0 if SE is well-calibrated).\n")
cat("  - Cov.95   = empirical coverage of the 95% CI (should be near 0.950).\n")
cat("  - Mean.SE / SE.Ratio / Cov.95 are reported as --- for methods without\n")
cat("    analytic standard errors (Joint_modeling_..._Liang, IIRR-weighting,\n")
cat("    IrregLong::Liang). Bootstrap inference would be needed for those.\n\n")


# ============================================================================
# 8. Save everything
# ============================================================================

dir.create(dirname(OUT_FILE), recursive = TRUE, showWarnings = FALSE)

saveRDS(
  list(
    sim_df      = sim_df,
    metrics     = metrics,
    settings = list(
      n_sim       = N_SIM,
      n_subj      = N_SUBJ,
      scenarios   = SELECTED_SCENARIOS,
      methods     = SELECTED_METHODS,
      true_gamma  = TRUE_GAMMA,
      true_alpha  = TRUE_ALPHA,
      true_beta   = TRUE_BETA,
      true_delta_template = TRUE_DELTA,
      true_theta_template = TRUE_THETA
    ),
    elapsed_min = elapsed
  ),
  file = OUT_FILE
)

cat(sprintf("Results saved to %s\n", OUT_FILE))

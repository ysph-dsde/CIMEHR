#' Add Realistic HbA1c Columns to Simulated EHR Data
#'
#' @description
#' Internal helper that takes the long-format \code{data.frame} produced by
#' \code{\link{sim_data_gen}} (whose outcome column \code{Y} is on a
#' standardized linear scale) and appends two extra columns: \code{HbA1c},
#' a clinically interpretable raw value in approximately the [4.5, 9] range,
#' and \code{log_HbA1c}, the natural-log transform used by analysis models.
#'
#' The mapping is
#' \deqn{\mathrm{log\_HbA1c} = \mathrm{log\_offset} +
#'       \mathrm{log\_scale} \cdot Y, \quad
#'       \mathrm{HbA1c} = \exp(\mathrm{log\_HbA1c}).}
#' With the default \code{log_offset = log(5.6)} and
#' \code{log_scale = 0.05}, a one-unit change in \code{Y} corresponds to
#' approximately a 5\% change in HbA1c, which is a clinically reasonable
#' effect size and keeps the population HbA1c distribution roughly in the
#' expected outpatient range.
#'
#' By default the original \code{Y} column is dropped, since the bundled
#' dataset documents and supports analysis on \code{log_HbA1c}. Set
#' \code{keep_Y = TRUE} to retain it.
#'
#' @param data A \code{data.frame} with a column named \code{Y}, typically
#'   \code{sim$long_data} from \code{\link{sim_data_gen}}.
#' @param log_offset Numeric; the constant added to \code{log_scale * Y} to
#'   produce \code{log_HbA1c}. Defaults to \code{log(5.6)} so that the
#'   population mean HbA1c sits near 5.6\% when \code{E[Y] = 0}.
#' @param log_scale Numeric; the multiplicative scaling on \code{Y} before
#'   exponentiation. Defaults to \code{0.05}.
#' @param keep_Y Logical; whether to retain the standardized \code{Y}
#'   column. Default \code{FALSE}.
#'
#' @return The input \code{data.frame} with new columns \code{HbA1c} (raw)
#'   and \code{log_HbA1c} (natural log of \code{HbA1c}). When the original
#'   \code{Y} is \code{NA}, both new columns are also \code{NA}.
#'
#' @keywords internal
#' @noRd
add_hba1c_columns <- function(data,
                              log_offset = log(5.6),
                              log_scale  = 0.05,
                              keep_Y     = FALSE) {
  if (!"Y" %in% names(data)) {
    stop("`data` must contain a column named `Y`.", call. = FALSE)
  }
  log_h <- log_offset + log_scale * data$Y
  data$log_HbA1c <- log_h
  data$HbA1c     <- exp(log_h)
  if (!keep_Y) data$Y <- NULL
  data
}

## ----  data-raw/sim_ehr_data.R  -----------------------------------------------
## Run once to (re)generate the example dataset shipped with the package.
##
## Design:
##   The bundled dataset is one replicate from the package's Monte Carlo
##   simulation framework (see simulations/simulation_study.R), specifically
##   the `continuous_phi0` scenario (focal Z = NSES, phi = 0). This keeps the
##   bundled dataset and the simulation study on a single, documented
##   data-generating spec.
##
##   - EHR-style covariate set: Age, Gender, Marital, Black, Hispanic, NSES.
##     Focal coefficient for the bundled dataset is on NSES (continuous).
##   - n = 2000, phi = 0 (no OU serial correlation).
##   - Time window [0, 120]; C is constant at time.end = 120 across subjects.
##   - Latent quantities (U, eta, q_*, b_*) are NOT included; they would not
##     be observable in real EHR data.
##   - Age is z-standardized (the `standardize = TRUE` covariate spec entry).
##     The raw-years `Age_raw` column that sim_data_gen() also produces is
##     dropped from the bundled data frame for a cleaner user-facing schema.
##   - The outcome is provided as a single column:
##       log_HbA1c - natural log of HbA1c, recommended for analysis models.
##                   Built from the standardized linear outcome Y via
##                   log_HbA1c = log(5.6) + 0.05 * Y.
##     The intermediate Y column and the raw HbA1c column are dropped from
##     the bundled data frame.
##
## True coefficients follow the simulation_study.R framework:
##   gamma (visiting): Age = 0.5 (X role); NSES = 0.3 (Z, focal);
##     background values for Gender, Marital, Black, Hispanic.
##   alpha (observation): Age = -0.3; NSES = -0.6; background for the rest.
##   beta (outcome, Y scale): Age = 0.5; NSES = -0.5; background for the
##     rest. The methods then operate on log_HbA1c, so the implied
##     log_HbA1c-scale betas are 0.05 * Y-scale betas.
##   Variance/scale parameters: sigma_zeta = 0.5, sigma_y = 0.7,
##     beta_t = 0.01, Sigma_b = diag(0.5, 0.8).
##   Seed: BASE_SEED + 1000000L + 1L = 21260509 (sim_id = 1 of the
##     continuous_phi0 scenario in simulation_study.R).
##
## Usage (from the package root):
##   source("data-raw/sim_ehr_data.R")
## ---------------------------------------------------------------------------

library(CIMEHR)

# -------------------------------------------------------------------------

ehr_covs <- list(
  Age      = list(type = "continuous", dist = "uniform",
                  min = 18, max = 99, standardize = TRUE),
  Gender   = list(type = "binary", prob = 0.5),
  Marital  = list(type = "binary", prob = 0.4),
  Black    = list(type = "binary", prob = 0.30),
  Hispanic = list(type = "binary", prob = 0.20),
  NSES     = list(type = "continuous", dist = "normal",
                  mean = 0, sd = 1, standardize = FALSE)
)

sim <- sim_data_gen(
  n          = 2000,
  time.start = 0,
  time.end   = 120,
  seed       = 21260509,
  covariates = ehr_covs,

  gamma = c(intercept = -2.8,
            Age       =  0.5,
            Gender    = -0.10,
            Marital   = -0.05,
            Black     =  0.15,
            Hispanic  = -0.05,
            NSES      =  0.3),
  sigma_zeta = 0.5,

  alpha = c(intercept =  0.2,
            Age       = -0.3,
            Gender    = -0.10,
            Marital   =  0.00,
            Black     =  0.05,
            Hispanic  = -0.05,
            NSES      = -0.6),
  obs_random_covs = "NSES",
  delta   = c(intercept = -0.3, NSES = -0.3),
  sigma_q = c(intercept =  0.5, NSES =  0.2),
  phi     = 0,

  beta = c(intercept = -2.0,
           Age       =  0.5,
           Gender    =  0.10,
           Marital   = -0.03,
           Black     =  0.15,
           Hispanic  =  0.05,
           NSES      = -0.5),
  beta_t              = 0.01,
  outcome_random_covs = "NSES",
  theta   = c(intercept = 0.4, NSES = 0.0),
  Sigma_b = matrix(c(0.5, 0, 0, 0.8), 2, 2),
  sigma_y = 0.7
)

# ---- Add HbA1c columns and finalize ----------------------------------------
# keep_Y = FALSE drops the intermediate standardized linear outcome from the
# bundled data frame; analyses use log_HbA1c (or HbA1c).
sim_ehr_data <- add_hba1c_columns(sim$long_data, keep_Y = FALSE)

# Drop the raw (un-standardized) Age column. The modeling functions use the
# standardized `Age` column; `Age_raw` was previously kept for descriptive
# summaries but is omitted from the bundled dataset for a cleaner schema.
sim_ehr_data$Age_raw <- NULL

# Drop the raw HbA1c column. Analyses use `log_HbA1c` (recommended) and the
# bundled-data schema documents that column only.
sim_ehr_data$HbA1c <- NULL

# Reorder columns
preferred_order <- c("id", "Age", "Gender", "Marital", "Black",
                     "Hispanic", "NSES", "time", "R",
                     "log_HbA1c", "C")
sim_ehr_data <- sim_ehr_data[, intersect(preferred_order, names(sim_ehr_data))]

# Attach true-parameter attribute AFTER column subsetting (which would otherwise
# strip user-added attributes from the data.frame).
attr(sim_ehr_data, "true_params") <- sim$params

# ---- Save ------------------------------------------------------------------
usethis::use_data(sim_ehr_data, overwrite = TRUE)

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
#' The original \code{Y} column is retained so that downstream models can
#' choose between fitting on \code{Y}, \code{log_HbA1c} (recommended for
#' parity with the manuscript) or \code{HbA1c}.
#'
#' @param data A \code{data.frame} with a column named \code{Y}, typically
#'   \code{sim$long_data} from \code{\link{sim_data_gen}}.
#' @param log_offset Numeric; the constant added to \code{log_scale * Y} to
#'   produce \code{log_HbA1c}. Defaults to \code{log(5.6)} so that the
#'   population mean HbA1c sits near 5.6\% when \code{E[Y] = 0}.
#' @param log_scale Numeric; the multiplicative scaling on \code{Y} before
#'   exponentiation. Defaults to \code{0.05}.
#' @param keep_Y Logical; whether to retain the standardized \code{Y}
#'   column. Default \code{TRUE}.
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
                              keep_Y     = TRUE) {
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
##   - EHR-style covariate set: Age, Gender, Marital, Black, Hispanic, NSES.
##     Focal coefficient for the bundled dataset is on NSES (continuous).
##   - n = 2000, phi = 0 (no OU serial correlation), seed = 2026.
##   - C is constant at time.end = 60 across subjects (no random censoring).
##   - Latent quantities (U, eta, q_*, b_*) are NOT included; they would not
##     be observable in real EHR data.
##   - Age has TWO columns:
##       Age      - z-standardized; this is what the modeling functions use
##       Age_raw  - original years (Uniform(18, 99) draw); for descriptive
##                  summaries only. Generated automatically by sim_data_gen()
##                  when a covariate is declared with `standardize = TRUE`.
##   - The outcome is provided in three columns:
##       Y         - standardized linear outcome from the DGP (analytic scale
##                   matching the manuscript's simulation Setting A.4)
##       HbA1c     - clinically interpretable raw value (~ 4.0-9.0)
##       log_HbA1c - natural log of HbA1c, recommended for analysis models
##
## True coefficient choices: beta_NSES = -0.25 mirrors the GIVEHR All of Us
## point estimate. Other betas are small-to-moderate, signed to match the
## cohort findings (Age +, Gender +, Black +, Hispanic small +, Marital small -).
## sigma_zeta = 1.0 and theta_NSES = 0.30 match the manuscript's Setting A.4
## ("full latent dependence") regime, where IP-only methods exhibit clear bias
## and CIMEHR recovers near-zero bias.
##
## Usage (from the package root):
##   source("data-raw/sim_ehr_data.R")
## ---------------------------------------------------------------------------

#library(CIMEHR)

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
  time.end   = 60,
  seed       = 2026,
  covariates = ehr_covs,

  gamma = c(intercept = -2.0,
            Age       =  0.10,
            Gender    = -0.10,
            Marital   = -0.05,
            Black     =  0.15,
            Hispanic  = -0.05,
            NSES      =  0.12),
  sigma_zeta = 1.0,

  alpha = c(intercept =  0.5,
            Age       =  0.10,
            Gender    = -0.10,
            Marital   =  0.00,
            Black     =  0.05,
            Hispanic  = -0.05,
            NSES      =  0.20),
  obs_random_covs = "NSES",
  delta   = c(intercept = -0.3, NSES = -0.5),
  sigma_q = c(intercept =  0.3, NSES =  0.2),
  phi     = 0,

  beta = c(intercept = -2.0,
           Age       =  0.05,
           Gender    =  0.10,
           Marital   = -0.03,
           Black     =  0.15,
           Hispanic  =  0.05,
           NSES      = -0.25),
  beta_t              = 0.01,
  outcome_random_covs = "NSES",
  theta   = c(intercept = 0.5, NSES = 0.30),
  Sigma_b = matrix(c(1, 0, 0, 2), 2, 2),
  sigma_y = 1
)

# ---- Add HbA1c columns and finalize ----------------------------------------
sim_ehr_data <- add_hba1c_columns(sim$long_data)
attr(sim_ehr_data, "true_params") <- sim$params

# Reorder columns so Age_raw sits next to Age for readability
preferred_order <- c("id", "Age", "Age_raw", "Gender", "Marital", "Black",
                     "Hispanic", "NSES", "time", "R", "Y", "HbA1c",
                     "log_HbA1c", "C")
sim_ehr_data <- sim_ehr_data[, intersect(preferred_order, names(sim_ehr_data))]

# ---- Save ------------------------------------------------------------------
usethis::use_data(sim_ehr_data, overwrite = TRUE)

## summarize_sim_ehr_data.R
## Purpose: produce a baseline / visit-pattern summary of the simulated EHR
## dataset bundled with the CIMEHR package.
##
## Usage:
##   source("summarize_sim_ehr_data.R")
##   summarize_sim_ehr_data(sim_ehr_data)
##
## If the package data are installed:
##   library(CIMEHR)
##   data(sim_ehr_data)
##   summarize_sim_ehr_data(sim_ehr_data)

summarize_sim_ehr_data <- function(dat,
                                   id_col       = "id",
                                   age_col      = "Age",
                                   age_raw_col  = "Age_raw",
                                   gender_col   = "Gender",
                                   marital_col  = "Marital",
                                   black_col    = "Black",
                                   hispanic_col = "Hispanic",
                                   nses_col     = "NSES",
                                   time_col     = "time",
                                   obs_col      = "R",
                                   outcome_col  = "HbA1c") {
  stopifnot(is.data.frame(dat))

  required <- c(id_col, gender_col, marital_col, black_col,
                hispanic_col, nses_col, time_col, obs_col)
  missing_required <- setdiff(required, names(dat))
  if (length(missing_required) > 0) {
    stop("The dataset is missing required columns: ",
         paste(missing_required, collapse = ", "), call. = FALSE)
  }

  pct <- function(x) 100 * mean(x == 1, na.rm = TRUE)

  iqr_text <- function(x) {
    q <- stats::quantile(x, probs = c(0.25, 0.50, 0.75),
                         na.rm = TRUE, names = FALSE)
    sprintf("%.1f [%.1f, %.1f]", q[2], q[1], q[3])
  }

  # ---- One row per subject for baseline covariates -------------------------
  baseline_cols <- intersect(c(id_col, age_col, age_raw_col, gender_col,
                               marital_col, black_col, hispanic_col, nses_col),
                             names(dat))
  baseline <- dat[!duplicated(dat[[id_col]]), baseline_cols, drop = FALSE]

  # ---- Visit and measurement counts per subject ----------------------------
  row_has_visit       <- !is.na(dat[[time_col]])
  row_has_measurement <- !is.na(dat[[obs_col]]) & dat[[obs_col]] == 1

  visits_by_id       <- tapply(row_has_visit,       dat[[id_col]], sum)
  measurements_by_id <- tapply(row_has_measurement, dat[[id_col]], sum)

  ids <- baseline[[id_col]]
  visits_by_id       <- visits_by_id[as.character(ids)]
  measurements_by_id <- measurements_by_id[as.character(ids)]
  visits_by_id[is.na(visits_by_id)]             <- 0
  measurements_by_id[is.na(measurements_by_id)] <- 0

  n_total            <- length(ids)
  n_with_visit       <- sum(visits_by_id >= 1)
  n_no_visit         <- sum(visits_by_id == 0)
  n_with_measurement <- sum(measurements_by_id >= 1)

  # ---- Age (years) ---------------------------------------------------------
  # Age in the modeling columns is z-standardized. The companion column
  # Age_raw stores the original years (Uniform(18, 99) draw) for descriptive
  # use only. If Age_raw is absent, fall back to a clear note.
  if (age_raw_col %in% names(baseline)) {
    age_years <- baseline[[age_raw_col]]
    age_mean_text  <- sprintf("%.1f", mean(age_years, na.rm = TRUE))
    age_sd_text    <- sprintf("%.1f", stats::sd(age_years, na.rm = TRUE))
    age_range_text <- sprintf("%.1f to %.1f",
                              min(age_years, na.rm = TRUE),
                              max(age_years, na.rm = TRUE))
  } else {
    age_mean_text  <- paste0("not available (only standardized ", age_col,
                             " in the dataset)")
    age_sd_text    <- "not available"
    age_range_text <- "not available"
  }

  # ---- White definition ----------------------------------------------------
  # Black and Hispanic are independent indicators in the simulator, so this
  # is a derived "non-Black, non-Hispanic" share rather than a true White %.
  nbnh_pct <- 100 * mean(baseline[[black_col]] == 0 &
                           baseline[[hispanic_col]] == 0, na.rm = TRUE)

  # ============================================================
  # Section 1: Sample size and baseline characteristics
  # ============================================================
  baseline_df <- data.frame(
    Characteristic = c(
      "Total patients",
      "Patients with >= 1 visit",
      "Patients with no visits",
      "Mean age (years)",
      "SD age (years)",
      "Age range (years)",
      "Male, %",
      "Married, %",
      "Black, %",
      "Hispanic, %",
      "Non-Black, non-Hispanic, %",
      "NSES mean",
      "NSES SD",
      "NSES range"
    ),
    Value = c(
      as.character(n_total),
      as.character(n_with_visit),
      as.character(n_no_visit),
      age_mean_text,
      age_sd_text,
      age_range_text,
      sprintf("%.1f", pct(baseline[[gender_col]])),
      sprintf("%.1f", pct(baseline[[marital_col]])),
      sprintf("%.1f", pct(baseline[[black_col]])),
      sprintf("%.1f", pct(baseline[[hispanic_col]])),
      sprintf("%.1f", nbnh_pct),
      sprintf("%.3f", mean(baseline[[nses_col]], na.rm = TRUE)),
      sprintf("%.3f", stats::sd(baseline[[nses_col]], na.rm = TRUE)),
      sprintf("%.3f to %.3f",
              min(baseline[[nses_col]], na.rm = TRUE),
              max(baseline[[nses_col]], na.rm = TRUE))
    ),
    stringsAsFactors = FALSE
  )

  # ============================================================
  # Section 2: HbA1c measurement and visit-pattern summary
  # ============================================================
  total_visits       <- sum(visits_by_id)
  total_measurements <- sum(measurements_by_id)
  meas_rate_pct      <- if (total_visits > 0) {
    100 * total_measurements / total_visits
  } else {
    NA_real_
  }

  measurement_df <- data.frame(
    Characteristic = c(
      "Patients with >= 1 HbA1c measurement",
      "Patients with >= 1 HbA1c measurement, %",
      "Visits per person, median [IQR]",
      "HbA1c measurements per person, median [IQR]",
      "Total visit rows",
      "Total observed HbA1c measurements",
      "Observed measurements among visits, %"
    ),
    Value = c(
      as.character(n_with_measurement),
      sprintf("%.1f", 100 * n_with_measurement / n_total),
      iqr_text(visits_by_id),
      iqr_text(measurements_by_id),
      as.character(total_visits),
      as.character(total_measurements),
      sprintf("%.1f", meas_rate_pct)
    ),
    stringsAsFactors = FALSE
  )

  # ---- Print ---------------------------------------------------------------
  cat("\nBaseline characteristics\n")
  cat(strrep("-", 60), "\n", sep = "")
  print(baseline_df, row.names = FALSE, right = FALSE)

  cat("\nVisit and HbA1c measurement patterns\n")
  cat(strrep("-", 60), "\n", sep = "")
  print(measurement_df, row.names = FALSE, right = FALSE)

  invisible(list(
    baseline           = baseline_df,
    measurement        = measurement_df,
    visits_by_person   = data.frame(
      id              = ids,
      n_visits        = as.integer(visits_by_id),
      n_measurements  = as.integer(measurements_by_id)
    )
  ))
}

# ---- Run when sourced ------------------------------------------------------
library(CIMEHR)
data(sim_ehr_data)
summarize_sim_ehr_data(sim_ehr_data)

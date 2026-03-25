## ----  data-raw/sim_ehr_data.R  -----------------------------------------------
## Run this script once to generate the example dataset shipped with the package.
## Place this file in data-raw/ and add data-raw to .Rbuildignore.
##
## Usage (from the package root):
##   source("data-raw/sim_ehr_data.R")
## ---------------------------------------------------------------------------

library(CIMEHR)

set.seed(2026)
sim <- sim_data_gen(n = 2000, seed = 2026)

# The long-format data frame used in examples and the vignette
sim_ehr_data <- sim$long_data

# Add the censoring-time column required by GIVEHR()
sim_ehr_data$C <- sim$params$time.end

# Store the true parameter values as an attribute so users can compare
attr(sim_ehr_data, "true_params") <- sim$params

usethis::use_data(sim_ehr_data, overwrite = TRUE)

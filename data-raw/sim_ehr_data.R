## ----  data-raw/sim_ehr_data.R  -----------------------------------------------
## Run this script once to generate the example dataset shipped with the package.
## Place this file in data-raw/ and add data-raw to .Rbuildignore.
##
## Parameter choices are documented below. The key design goal is a dataset
## where the two informative-presence bias channels are clearly visible:
##   (1) Informative visiting: subjects visit more/less based on the shared
##       latent factor U, which also drives the outcome.
##   (2) Informative observation: whether the outcome is recorded at a visit
##       also depends on U via the negative delta link.
## Naive methods (LMM, SummaryStat) ignore both channels and are noticeably
## biased; GIVEHR recovers the truth.
##
## Usage (from the package root):
##   source("data-raw/sim_ehr_data.R")
## ---------------------------------------------------------------------------

library(CIMEHR)

## ---- Data-generating parameters -------------------------------------------
##
## Visiting process: lambda_i(t) = eta_i * exp(gamma' X^V)
##   eta_i = exp(mu0 + sigma_zeta * U_i),  U_i ~ N(0,1)
##
##   gamma_intercept = -2.0  (slightly raised from default -2.2; ~13 visits/subject)
##   gamma_Z = 1.0, gamma_X = 1.0  (unchanged: strong covariate effect on visit rate)
##   sigma_zeta = 1.2  (raised from 1.0: stronger frailty variation amplifies the
##                       correlation between visit frequency and the latent factor U,
##                       which is the primary driver of informative-visiting bias)
##
## Observation process: R_ij = I(alpha'X + q_i'Z^O + xi(t) > 0)
##   q_i | U_i ~ N(delta * U_i, Sigma_q),  xi(t) ~ OU(phi=0.1)
##
##   alpha_intercept = 0.5  (lowered from 1.5: reduces marginal observation rate
##                            from ~92% to ~67%, making missingness meaningful)
##   alpha_Z = -0.5, alpha_X = -0.3  (unchanged)
##   delta_intercept = -0.5  (strengthened from -0.3: high-U subjects visit more
##   delta_Z = -0.7           but are less likely to have Y recorded, creating
##                             informative observation bias on top of visiting bias)
##   sigma_q = (0.3, 0.2), phi = 0.1  (unchanged)
##
## Outcome process: Y_ij = beta'X + beta_t*t + b_i'Z^Y + epsilon
##   b_i | U_i ~ N(theta * U_i, Sigma_b)
##
##   beta = (-2, Z=-0.5, X=0.5)  (unchanged: these are the target estimands)
##   theta_intercept = 0.5  (unchanged)
##   theta_Z = 0.5  (raised from 0.3: stronger link between U and outcome,
##                    amplifying the omitted-variable bias in naive methods)
##   Sigma_b = diag(1, 2), sigma_y = 1  (unchanged)
##
## Net effect on naive-method bias (approximate):
##   The proxy bias for beta_Z scales with theta_Z * Corr(visits, U).
##   Default:  0.3 * 0.667 = 0.200
##   Proposed: 0.5 * 0.732 = 0.366  (~1.8x amplification)
##
## Dataset properties (approximate):
##   Expected visits per subject:     13
##   Expected observation rate:       67%
##   Expected observed Y per subject:  9

sim <- sim_data_gen(
  n    = 2000,
  seed = 101,

  # ---- Visiting process ----
  gamma      = c(intercept = -2.0, Z = 1.0, X = 1.0),
  sigma_zeta = 1.2,

  # ---- Observation process ----
  alpha   = c(intercept =  0.5, Z = -0.5, X = -0.3),
  delta   = c(intercept = -0.5, Z = -0.7),
  sigma_q = c(intercept =  0.3, Z =  0.2),
  phi     = 0.1,

  # ---- Outcome process ----
  beta    = c(intercept = -2, Z = -0.5, X = 0.5),
  beta_t  = 0.1,
  theta   = c(intercept = 0.5, Z = 0.5),
  Sigma_b = matrix(c(1, 0, 0, 2), 2, 2),
  sigma_y = 1
)

# The long-format data frame used in examples and the vignette
sim_ehr_data      <- sim$long_data
sim_ehr_data$C    <- sim$params$time.end
attr(sim_ehr_data, "true_params") <- sim$params

usethis::use_data(sim_ehr_data, overwrite = TRUE)

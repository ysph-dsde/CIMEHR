CIMEHR
================

<!-- badges: start -->

[![R-CMD-check](https://github.com/ysph-dsde/CIMEHR/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/ysph-dsde/CIMEHR/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

> A unified semiparametric joint modeling framework that simultaneously
> characterizes the visiting, biomarker observation, and longitudinal
> outcome processes.

# Overview

**CIMEHR**, which stands for **C**linical **I**nformative **M**issingness for **E**lectronic 
**H**ealth **R**ecord data, fits semiparametric joint models for longitudinal EHR 
data that simultaneously correct for two sources of informative presence bias: 
an **informative visiting process** (when patients visit) and an 
**informative observation process** (whether outcomes are recorded at visits).

The primary method is **`GIVEHR()`** (Generalized Informative Visiting and 
observation processes in Electronic Health Records), a three-stage estimator that:

1.  models *when* patients visit via partial likelihood with a
    log-normal frailty (visiting process),
2.  models *whether* the outcome is recorded at a visit via probit
    regression with shared-latent-factor-linked random effects
    (observation process),
3.  estimates the longitudinal outcome coefficients via weighted least
    squares with risk-set centering (outcome process).

All three stages are connected through a shared latent factor. Sandwich
variance estimators provide robust inference throughout.

## Installation

``` r
# install.packages("devtools")
devtools::install_github("ysph-dsde/CIMEHR")
```

## Quick example

The package ships with a pre-simulated dataset `sim_ehr_data` (2000
subjects) so you can get started immediately:

``` r
library(CIMEHR)

# Use the bundled dataset (already includes censoring time C)
data(sim_ehr_data)

# Fit the three-stage joint model
fit <- GIVEHR(data = sim_ehr_data)

# Outcome fixed effects with SE and 95% CI
print(fit)

# Full summary of all three processes
summary(fit)
```

To simulate fresh data with custom parameters, use `sim_data_gen()`:

``` r
dat <- sim_data_gen(n = 500, seed = 42)

# Add censoring time (not included by sim_data_gen)
dat$long_data$C <- 60

fit <- GIVEHR(data = dat$long_data)
```

### Extracting results programmatically

``` r
# Outcome coefficients (beta) - default
coef(fit)
coef(fit, type = "visiting")
coef(fit, type = "observation")

# Variance-covariance matrix and confidence intervals
vcov(fit)
confint(fit)

# Direct access to components
fit$longitudinal_outcome$beta_Y
fit$longitudinal_outcome$se_beta_naive
fit$visiting_process$gamma_V
fit$observation_process$alpha_fixed
```

## GIVEHR variants

Two extensions handle temporal dependence in the observation process:

- **`GIVEHR_GH()`**: Gauss-Hermite quadrature with online Laplace
  filtering. Supports time-varying observation dependence via a lagged
  indicator.
- **`GIVEHR_OU_PCL()`**: Ornstein-Uhlenbeck pairwise composite
  likelihood. Captures serial correlation in the observation process.

``` r
# Gauss-Hermite variant
fit_gh <- GIVEHR_GH(data = sim_ehr_data, gh_points_U = 7, gh_points_q = 5)

# OU + Pairwise composite likelihood variant
fit_ou <- GIVEHR_OU_PCL(data = sim_ehr_data, pair_type = "adjacent")
```

## Benchmark methods

The package includes several comparison methods for simulation and
benchmark studies:

| Function | Approach |
|----|----|
| `EHRJoint()` | Joint model without separate observation process |
| `JMVL_Liang()` | Liang–Lu–Ying joint estimating equations |
| `IIRR_weighting()` | Inverse intensity rate ratio (IIRR) weighting |
| `IIRR_balancing()` | Inverse intensity rate ratio (IIRR) with balancing weights |
| `LMM()` | Linear mixed model (standard, observed-adjusted, visit-adjusted variants) |
| `LI_IP()` | Linear increment imputation + informative visit (IP) adjustment |
| `MI_IP()` | Multiple imputation + informative visit (IP) adjustment |
| `PairwiseLikelihood()` | Pairwise composite likelihood ranking |
| `SummaryStat()` | Summary-statistic regression (naive baseline) |

## Getting help

- Run `vignette("getting-started", package = "CIMEHR")` for a complete
  tutorial.
- Open an issue on GitHub for bugs or feature requests.

## References

Yang, C.-H., Shi, X., and Mukherjee, B. (2026). “Joint Modeling of
Longitudinal EHR Data with Shared Random Effects for Informative
Visiting and Observation Processes.” *arXiv:2602.15374*.

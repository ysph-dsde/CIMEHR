# SIMMER

<!-- badges: start -->
[![R-CMD-check](https://github.com/ysph-dsde/SIMMER/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/ysph-dsde/SIMMER/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

> A unified semiparametric joint modeling framework that
> simultaneously characterizes the visiting, biomarker observation,
> and longitudinal outcome processes.

---
## Overview

**SIMMER** fits semiparametric joint models for longitudinal EHR data
that simultaneously correct for two sources of informative presence
bias: an **informative visiting process** (when patients visit) and an
**informative observation process** (whether outcomes are recorded at
visits).

The primary method is **`GIVEHR()`** (Generalized Informative Visiting
and Examination in EHR Data), a three-stage estimator that:

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
devtools::install_github("ysph-dsde/SIMMER")
```

## Quick example

``` r
library(SIMMER)

# Simulate EHR data with informative visiting and observation
dat <- sim_data_gen(n = 500, seed = 42)

# Add censoring time (required by GIVEHR)
dat$long_data$C <- 60

# Fit the three-stage joint model
fit <- GIVEHR(data = dat$long_data)

# Outcome fixed effects with SE and 95% CI
fit

# Full summary of all three processes
summary(fit)
```

### Extracting results programmatically

``` r
# Outcome coefficients (beta)
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

- **`GIVEHR_GH()`**: Gauss–Hermite quadrature with online Laplace
  filtering. Supports time-varying observation dependence via a lagged
  indicator.
- **`GIVEHR_OU_PCL()`**: Ornstein–Uhlenbeck pairwise composite
  likelihood. Captures serial correlation in the observation process.

``` r
# Gauss-Hermite variant
fit_gh <- GIVEHR_GH(data = dat$long_data, gh_points_U = 7, gh_points_q = 5)

# OU + Pairwise composite likelihood variant
fit_ou <- GIVEHR_OU_PCL(data = dat$long_data, pair_type = "adjacent")
```

## Benchmark methods

The package includes several comparison methods for simulation and benchmark studies:

| Function | Approach |
|---|---|
| `EHRJoint()` | Joint model without separate observation process |
| `JMVL_IP()` | Liang--Lu--Ying joint estimating equations |
| `IIRR_weighting()` | Inverse intensity rate ratio (IIRR) weighting |
| `IIRR_balancing()` | Inverse intensity rate ratio (IIRR) with balancing weights |
| `LMM()` | Linear mixed model (standard, observed-adjusted, visit-adjusted variants) |
| `LI_IP()` | Linear increment imputation + informative visit (IP) adjustment |
| `MI_IP()` | Multiple imputation + informative visit (IP) adjustment |
| `PairwiseLikelihood()` | Pairwise composite likelihood ranking |
| `SummaryStat()` | Summary-statistic regression (naive baseline) |

## Getting help

- Run `vignette("getting-started", package = "SIMMER")` for a complete
  tutorial.
- Open an issue on GitHub for bugs or feature requests.

## References

Yang, C.-H., Shi, X., and Mukherjee, B. (2026). "Joint Modeling of
Longitudinal EHR Data with Shared Random Effects for Informative
Visiting and Observation Processes." *arXiv:2602.15374*.

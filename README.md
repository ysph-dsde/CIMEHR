
<!-- badges: start -->
[![R-CMD-check](https://github.com/ysph-dsde/CIMEHR/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/ysph-dsde/CIMEHR/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

> Semiparametric and benchmark methods for longitudinal EHR data with informative visiting and informative observation processes.

# Overview

**CIMEHR** (**C**linical **I**nformative **M**issingness for **E**lectronic **H**ealth **R**ecord data) provides methods for longitudinal EHR analyses where patients visit irregularly and outcomes may be selectively measured at visits. The primary methods jointly characterize three linked processes:

1. the **visiting process**: when a patient has an EHR visit;
2. the **observation process**: whether the outcome/biomarker is recorded at a visit;
3. the **longitudinal outcome process**: the outcome trajectory among the target population.

The main estimator is `CIMEHR()`. Two extensions, `CIMEHR_timevarying_ghq()` and `CIMEHR_timevarying_ou()`, allow more flexible time-varying/serially correlated observation-process structures.

## Installation

``` r
# install.packages("devtools")
devtools::install_github("ysph-dsde/CIMEHR")
```

## Data format

Most functions in **CIMEHR** expect **long-format longitudinal data**, with one row per subject-visit or subject-time record. The primary CIMEHR estimators require at least these columns:

| Column role | Default name | Meaning |
|:--|:--|:--|
| Subject identifier | `id` | Unique patient/subject ID |
| Visit time | `time` | Visit time; `NA` can be used for subjects/rows without a visit |
| Outcome | `Y` | Longitudinal outcome; usually `NA` when not observed |
| Observation indicator | `R` | 1 if `Y` was observed at that visit, 0 otherwise |
| Censoring/follow-up end | `C` | End of follow-up for each subject |
| Covariates | e.g., `Z`, `X` | Baseline or model covariates |

The package methods are designed for long/panel data. Some benchmark methods use only observed rows, but the input should still be organized as long data so visit timing, observation indicators, and outcome measurements are aligned.

### Converting short/wide data to long data

If your data are stored with one row per subject and repeated outcome/time columns, convert them to long format before fitting a model:

``` r
library(dplyr)
library(tidyr)

wide_dat <- tibble::tibble(
  id = 1:2,
  Z = c(0, 1),
  X = c(0.2, -0.4),
  C = c(60, 60),
  time_1 = c(5, 8),  time_2 = c(20, NA), time_3 = c(45, 50),
  Y_1 = c(1.2, NA), Y_2 = c(1.6, NA), Y_3 = c(2.0, 1.1)
)

long_dat <- wide_dat %>%
  pivot_longer(
    cols = matches("^(time|Y)_\\\\d+$"),
    names_to = c(".value", "visit"),
    names_pattern = "(time|Y)_(\\\\d+)"
  ) %>%
  mutate(R = as.integer(!is.na(Y))) %>%
  arrange(id, time)
```

## Quick example

The package includes a simulated long-format dataset, `sim_ehr_data`.

``` r
library(CIMEHR)

data("sim_ehr_data")

fit <- CIMEHR(data = sim_ehr_data)
print(fit)
summary(fit)
```

To simulate new data:

``` r
dat <- sim_data_gen(n = 500, seed = 42)
fit <- CIMEHR(data = dat$long_data)
```

## Method comparison


<table>
  <thead>
    <tr>
      <th align="left">Method</th>
      <th>IP</th>
      <th>IO</th>
      <th align="left">Shared Latent</th>
      <th align="left">Adjustment Mechanism</th>
      <th align="left">Package</th>
    </tr>
  </thead>
  <tbody>
    <tr><td colspan="6"><em>Outcome-only</em></td></tr>
    <tr>
      <td>Summary statistics</td>
      <td>✗</td>
      <td>✗</td>
      <td>None</td>
      <td>Reduces the longitudinal series to a single scalar (e.g., min, mean, median, max) for regression.</td>
      <td>—</td>
    </tr>
    <tr>
      <td>Standard LMM<br><sub>(Laird and Ware, 1982)</sub></td>
      <td>✗</td>
      <td>✗</td>
      <td>None</td>
      <td>Conditions on observed covariates (valid only under MAR).</td>
      <td>lme4</td>
    </tr>
    <tr>
      <td>VA-LMM<br><sub>(Goldstein et al., 2016)</sub></td>
      <td>✓</td>
      <td>✗</td>
      <td>None</td>
      <td>Includes <em>visit count</em> as a fixed covariate proxy.</td>
      <td>lme4</td>
    </tr>
    <tr>
      <td>OA-LMM</td>
      <td>✗</td>
      <td>✓</td>
      <td>None</td>
      <td>Includes <em>observation count</em> as a fixed covariate proxy.</td>
      <td>lme4</td>
    </tr>

    <tr><td colspan="6"><em>IP-only</em></td></tr>
    <tr>
      <td>Liang<br><sub>(Liang et al., 2009)</sub></td>
      <td>✓</td>
      <td>✗</td>
      <td>Visiting &amp; Outcome</td>
      <td>Joint likelihood estimation via shared frailty linking visiting and outcome.</td>
      <td>IrregLong</td>
    </tr>
    <tr>
      <td>JMVL-G<br><sub>(Gasparini et al., 2020)</sub></td>
      <td>✓</td>
      <td>✗</td>
      <td>Visiting &amp; Outcome</td>
      <td>Joint likelihood modeling inter-visit times via shared random effects.</td>
      <td>merlin</td>
    </tr>
    <tr>
      <td>IIRR-Weighting<br><sub>(Burvzrova and Lumley, 2007)</sub></td>
      <td>✓</td>
      <td>✗</td>
      <td>None</td>
      <td>Inverse intensity weighting based on visiting intensity.</td>
      <td>—</td>
    </tr>
    <tr>
      <td>IIRR-Balanced<br><sub>(Yiu and Su, 2025)</sub></td>
      <td>✓</td>
      <td>✗</td>
      <td>None</td>
      <td>Inverse intensity weighting adjusted for covariate balance.</td>
      <td>—</td>
    </tr>
    <tr>
      <td>Pairwise Likelihood (PairLik)<br><sub>(Chen et al., 2015)</sub></td>
      <td>✓</td>
      <td>✗</td>
      <td>None</td>
      <td>Constructs a likelihood from all within-subject pairs of observations; conditioning on observed visit times eliminates the visit intensity from the likelihood.</td>
      <td>—</td>
    </tr>

    <tr><td colspan="6"><em>Imputation + IP-only</em></td></tr>
    <tr>
      <td>Multiple Imputation (MI)<br><sub>(Rubin, 1987)</sub></td>
      <td>✓</td>
      <td>✓</td>
      <td>Determined by IP-only</td>
      <td>Multiple imputation based on posterior draws, followed by IP-only methods.</td>
      <td>mice</td>
    </tr>
    <tr>
      <td>Linear Increment Methods (LI)<br><sub>(Aalen and Gunnes, 2010)</sub></td>
      <td>✓</td>
      <td>✓</td>
      <td>Determined by IP-only</td>
      <td>Deterministic linear interpolation between observed visits to impute missing values, followed by IP-only methods.</td>
      <td>slim</td>
    </tr>

    <tr><td colspan="6"><em>IP+IO</em></td></tr>
    <tr>
      <td>Adapted-Liang<br><sub>(Liang et al., 2009)</sub></td>
      <td>✓</td>
      <td>✓</td>
      <td>Visiting &amp; Outcome</td>
      <td>Joint model structure defined, but observation coefficients are constrained to zero.</td>
      <td>—</td>
    </tr>
    <tr>
      <td>EHRJoint<br><sub>(Du et al., 2025)</sub></td>
      <td>✓</td>
      <td>✓</td>
      <td>Visiting &amp; Outcome</td>
      <td>Tripartite model; observation process modeled via covariates only (no latent link).</td>
      <td>—</td>
    </tr>
    <tr>
      <td>Proposed<br><sub>(this package)</sub></td>
      <td>✓</td>
      <td>✓</td>
      <td>All three processes</td>
      <td>Tripartite model; full shared latent structure links Visiting, Observation, and Outcome.</td>
      <td>CIMEHR</td>
    </tr>
  </tbody>
</table>

IP = informative visiting / informative presence adjustment; IO = informative observation adjustment.
The “Proposed” row corresponds to <code>CIMEHR()</code>; <code>CIMEHR_timevarying_ghq()</code> and <code>CIMEHR_timevarying_ou()</code> extend this framework with more flexible observation-process modeling.

## Primary CIMEHR methods

``` r
fit_base <- CIMEHR(
  data = sim_ehr_data,
  covars_visit_XV = c("Z", "X"),
  covars_outcome_fixed_XY = c("Z", "X"),
  covars_outcome_random_link_ZY = "Z",
  covars_obs_fixed_XO = c("Z", "X"),
  covars_obs_random_link_ZO = "Z"
)

fit_ghq <- CIMEHR_timevarying_ghq(data = sim_ehr_data,
                                  gh_points_U = 7,
                                  gh_points_q = 5)

fit_ou <- CIMEHR_timevarying_ou(data = sim_ehr_data,
                                pair_type = "adjacent")
```

The `summary()` output for CIMEHR-family fits reports the frailty model and the observation-process link. For example, `CIMEHR()` and `CIMEHR_timevarying_ghq()` use a probit observation model, while `CIMEHR_timevarying_ou()` uses a probit observation model with OU serial correlation estimated by pairwise composite likelihood.

## Extracting results

``` r
coef(fit_base)
coef(fit_base, type = "visiting")
coef(fit_base, type = "observation")

vcov(fit_base)
confint(fit_base)

fit_base$observation_process$delta_latent_random_link
fit_base$longitudinal_outcome$theta_Y
```

## Bootstrap inference

Use `bootstrap()` for methods without analytic standard errors or as a robustness check for CIMEHR-family methods:

``` r
fit_iirr <- Inverse_intensity_rate_ratio_weighting(sim_ehr_data,
                                                   visit_covs = c("Z", "X"))
bs <- bootstrap(fit_iirr, data = sim_ehr_data, B = 200, seed = 1)
print(bs)
summary(bs)
```

Model-specific arguments can be passed through `covars`:

``` r
bs_cimehr <- bootstrap(
  fit_base,
  data = sim_ehr_data,
  B = 200,
  seed = 1,
  covars = list(
    covars_visit_XV = c("Z", "X"),
    covars_outcome_fixed_XY = c("Z", "X"),
    covars_outcome_random_link_ZY = "Z",
    covars_obs_fixed_XO = c("Z", "X"),
    covars_obs_random_link_ZO = "Z"
  )
)
```

## Getting help

Run `vignette("getting-started", package = "CIMEHR")` for a more detailed tutorial.

## References

Yang, C.-H., Shi, X., and Mukherjee, B. (2026). "Joint Modeling of Longitudinal EHR Data with Shared Random Effects for Informative Visiting and Observation Processes." *arXiv:2602.15374*.

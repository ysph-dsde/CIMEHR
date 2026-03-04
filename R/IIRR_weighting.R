#' Inverse Intensity Rate Ratio Weighting Estimator
#'
#' @description
#' Estimates marginal regression coefficients for a longitudinal outcome by
#' weighting each observation by the inverse of its estimated visit
#' intensity. The visiting process is modelled via a Cox-type partial
#' likelihood score equation, and the outcome coefficients are obtained by
#' solving inverse-intensity-weighted estimating equations.
#'
#' Only rows with observed outcomes (\code{R == 1} or \code{R} is
#' \code{NA}) are used in both the visiting-process and outcome
#' estimation.
#'
#' @param data A \code{data.frame} in long (panel) format with one row per
#'   subject-time combination. Must contain columns for subject identifier,
#'   visit time, outcome, observation indicator, and all specified
#'   covariates.
#' @param id Character string giving the name of the subject identifier
#'   column. Default is \code{"id"}.
#' @param time Character string giving the name of the visit-time column.
#'   Default is \code{"time"}.
#' @param outcome Character string giving the name of the outcome column.
#'   Default is \code{"Y"}.
#' @param obs_indicator Character string giving the name of the binary
#'   observation indicator column (\code{1} = outcome observed at visit).
#'   Default is \code{"R"}.
#' @param visit_covs Character vector of covariate column names for the
#'   visiting-process intensity model. These enter the exponent of the
#'   multiplicative intensity. Default is \code{c("Z", "X")}.
#' @param formula A one-sided \code{\link[stats]{formula}} specifying the
#'   design matrix for the longitudinal outcome model (e.g.,
#'   \code{~ Z + X + time}). An intercept is included by default. Default
#'   is \code{~ Z + X + time}.
#'
#' @return A named list with two components:
#'
#' \describe{
#'   \item{\code{gamma_hat}}{Named numeric vector of visiting-process
#'     regression coefficients, estimated by solving the partial-likelihood
#'     score equation with \code{\link[nleqslv]{nleqslv}}.}
#'   \item{\code{beta_hat}}{Named numeric vector of outcome model
#'     coefficients, estimated from inverse-intensity-weighted estimating
#'     equations solved with \code{\link[nleqslv]{nleqslv}}.}
#' }
#'
#' @details
#' The estimator proceeds in three steps:
#'
#' \strong{Step 1.} The visiting-process coefficients \eqn{\gamma} are
#' estimated by solving the partial-likelihood score equation
#' \deqn{\sum_j \bigl[Z_{ij} - \bar Z(\gamma; t_j)\bigr] = 0,}
#' where \eqn{\bar Z(\gamma; t) = \sum_i Z_i \exp(\gamma^\top Z_i) /
#' \sum_i \exp(\gamma^\top Z_i)} is the risk-set mean.
#'
#' \strong{Step 2.} Inverse intensity weights are computed as
#' \eqn{w_i = \exp(\hat\gamma^\top Z_i)}.
#'
#' \strong{Step 3.} The outcome coefficients \eqn{\beta} are estimated by
#' solving the weighted estimating equation
#' \deqn{\sum_{i,j} w_i^{-1} X_{ij} (Y_{ij} - X_{ij}^\top \beta) = 0.}
#'
#' @note
#' This function requires the \pkg{nleqslv} package. It will attempt to
#' install the package automatically if it is not available.
#'
#' No standard errors are returned. For inference, consider using
#' \code{\link{IIRR_balancing}} which provides balancing-weights
#' estimates, or bootstrap the entire procedure.
#'
#' @references
#' Petra Bůžková and Thomas Lumley (2007). Longitudinal data analysis
#' for generalized linear models with follow-up dependent on outcome-related variables.
#' \emph{The Canadian Journal of Statistics}, 35(4), 485--500. \doi{10.1002/cjs.5550350402}.
#'
#' @examples
#' \dontrun{
#' dat <- sim_data_gen(n = 500, seed = 42)
#'
#' fit <- IIRR_weighting(
#'   data      = dat$long_data,
#'   visit_covs = c("Z", "X"),
#'   formula   = ~ Z + X + time
#' )
#' fit$gamma_hat
#' fit$beta_hat
#' }
#'
#' @importFrom stats model.frame model.matrix na.pass
#' @export
IIRR_weighting <- function(data,
                           id = "id",
                           time = "time",
                           outcome = "Y",
                           obs_indicator = "R",
                           visit_covs = c("Z", "X"),
                           formula = ~ Z + X + time) {

  if (!requireNamespace("nleqslv", quietly = TRUE)) {
    stop("Package 'nleqslv' is required. Install it with: install.packages('nleqslv')",
         call. = FALSE)
  }

  # Filter to observed data only (R == 1 or R is NA)
  R_vec <- data[[obs_indicator]]
  data <- data[R_vec == 1 | is.na(R_vec), ]

  # Visiting process variables
  V <- as.matrix(data[, c(id, visit_covs)])
  VV <- as.matrix(V[!duplicated(V), -1])
  V <- as.matrix(data[, visit_covs])

  y <- data[[outcome]]

  # ============================================
  # Step 1: Estimate gamma (visiting process)
  # ============================================
  gamma_function <- function(gamma) {
    gamma <- matrix(gamma, ncol = 1)
    V_denom <- sum(exp(VV %*% gamma))
    V_numer <- colSums(VV * as.numeric(exp(VV %*% gamma)))
    obj <- colSums(V) - (V_numer / V_denom) * nrow(V)
    return(obj)
  }

  gamma_init <- rep(1, ncol(V))
  gamma_solve <- nleqslv::nleqslv(x = gamma_init, fn = gamma_function, control = list(trace = 0))
  gamma_hat <- gamma_solve$x
  names(gamma_hat) <- visit_covs

  # ============================================
  # Step 2: Compute weights
  # ============================================
  weights <- as.vector(exp(V %*% gamma_hat))

  # ============================================
  # Step 3: Build design matrix from formula
  # ============================================
  mf <- model.frame(formula, data = data, na.action = na.pass)
  X_GEE <- model.matrix(formula, data = mf)

  # ============================================
  # Step 4: Estimate beta (weighted estimating equation)
  # ============================================
  Y_function <- function(beta) {
    res <- as.vector((1 / weights) * (y - X_GEE %*% beta))
    obj <- colSums(na.omit(X_GEE * res))
    return(obj)
  }

  beta_init <- rep(0.1, ncol(X_GEE))
  beta_solve <- nleqslv::nleqslv(x = beta_init, fn = Y_function, control = list(trace = 0))
  beta_hat <- beta_solve$x
  names(beta_hat) <- colnames(X_GEE)

  return(list(
    gamma_hat = gamma_hat,
    beta_hat = beta_hat
  ))
}

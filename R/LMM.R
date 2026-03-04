#' Fit Linear Mixed Model for Longitudinal Data
#'
#' @param data A data frame containing longitudinal data.
#' @param id Character string for subject ID column. Default is "id".
#' @param fixed Formula for fixed effects (e.g., Y ~ Z + X + time).
#' @param random Formula for random effects (e.g., ~ 1 + Z | id).
#' @param obs_indicator Column name for observation indicator (keep rows where value != 0).
#'        Default is NULL (use all rows).
#' @param visit_adjust Type of visit count adjustment: NULL (none), "OA" (observed-adjusted,
#'        count after filtering), or "VA" (visit-adjusted, count before filtering).
#'        When specified, n_visits is automatically added to the fixed effects formula.
#' @param method Estimation method: "REML" or "ML". Default is "REML".
#' @param optimizer Optimizer: "nlminb" or "optim". Default is "nlminb".
#'
#' @return Coefficient matrix from summary(lme(...)).
#'
#' @examples
#' # Standard LMM
#' LMM(data, fixed = Y ~ Z + X + time, random = ~ 1 | id, obs_indicator = "R")
#'
#' # OA-LMM (n_visits automatically added to formula)
#' LMM(data, fixed = Y ~ Z + X + time, random = ~ 1 | id,
#'         obs_indicator = "R", visit_adjust = "OA")
#'
#' # VA-LMM (n_visits automatically added to formula)
#' LMM(data, fixed = Y ~ Z + X + time, random = ~ 1 | id,
#'         obs_indicator = "R", visit_adjust = "VA")
#'
#' @export
LMM <- function(data,
                    id = "id",
                    fixed = Y ~ 1,
                    random = ~ 1 | id,
                    obs_indicator = NULL,
                    visit_adjust = NULL,
                    method = "REML",
                    optimizer = "nlminb") {

  work_data <- data

  # VA: add visit count before filtering
  if (!is.null(visit_adjust) && toupper(visit_adjust) == "VA") {
    work_data$n_visits <- ave(seq_len(nrow(work_data)), work_data[[id]], FUN = seq_along)
  }

  # Filter by observation indicator
  if (!is.null(obs_indicator)) {
    work_data <- work_data[work_data[[obs_indicator]] != 0, ]
  }

  # OA: add visit count after filtering
  if (!is.null(visit_adjust) && toupper(visit_adjust) == "OA") {
    work_data$n_visits <- ave(seq_len(nrow(work_data)), work_data[[id]], FUN = seq_along)
  }

  # Add n_visits to formula if visit_adjust is specified
  if (!is.null(visit_adjust)) {
    fixed <- update(fixed, . ~ . + n_visits)
  }

  # Remove NA
  work_data <- na.omit(work_data)

  # Control options
  ctrl <- nlme::lmeControl(opt = optimizer, msMaxIter = 5000)
  if (optimizer == "optim") ctrl$msMaxEval <- 5000

  # Fit
  fit <- nlme::lme(fixed, data = work_data, random = random, method = method, control = ctrl)

  return(summary(fit)$tTable)
}

#tmp = LMM(dat$long_data, fixed = Y ~ Z + X + time, random = ~ 1 + Z| id,
#        obs_indicator = "R", visit_adjust = "VA")
#tmp

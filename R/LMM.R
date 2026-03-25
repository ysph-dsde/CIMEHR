#' Fit a Linear Mixed Model (LMM) for Longitudinal Data
#'
#' @description
#' Fits a linear mixed-effects model (LMM) to longitudinal data using
#' \code{\link[nlme]{lme}}, with optional filtering by an observation
#' indicator and optional visit-count adjustment. This provides a
#' convenient wrapper for standard LMM, Observation-Adjusted LMM (OA-LMM),
#' and Visit-Adjusted LMM (VA-LMM) estimators commonly used as
#' benchmarks in informative-visiting analyses.
#'
#' @param data A \code{data.frame} in long (panel) format containing
#'   longitudinal data with one row per subject-time combination.
#' @param id Character string giving the name of the subject identifier
#'   column. This is used for the visit-count adjustment (if requested)
#'   and is referenced in the random-effects formula. Default is
#'   \code{"id"}.
#' @param fixed A \code{\link[stats]{formula}} specifying the fixed-effects
#'   structure (e.g., \code{Y ~ Z + X + time}). When \code{visit_adjust}
#'   is non-\code{NULL}, the term \code{n_visits} is appended
#'   automatically. Default is \code{Y ~ 1}.
#' @param random A one-sided \code{\link[stats]{formula}} specifying the
#'   random-effects structure using \pkg{nlme} syntax (e.g.,
#'   \code{~ 1 + Z | id}). Default is \code{~ 1 | id}.
#' @param obs_indicator Character string giving the name of the
#'   observation indicator column. Rows where this column equals \code{0}
#'   are dropped before fitting. Set to \code{NULL} to use all rows.
#'   Default is \code{NULL}.
#' @param visit_adjust Character string specifying the type of visit-count
#'   adjustment. One of:
#'   \describe{
#'     \item{\code{NULL}}{(Default.) No adjustment; a standard LMM is
#'       fitted.}
#'     \item{\code{"OA"}}{Observed-adjusted: the cumulative visit count
#'       \code{n_visits} is computed \emph{after} filtering by
#'       \code{obs_indicator}. This adjusts for the number of observed
#'       (retained) visits.}
#'     \item{\code{"VA"}}{Visit-adjusted: the cumulative visit count
#'       \code{n_visits} is computed \emph{before} filtering by
#'       \code{obs_indicator}. This adjusts for the total number of
#'       visits regardless of whether the outcome was observed.}
#'   }
#'   When either \code{"OA"} or \code{"VA"} is specified, the variable
#'   \code{n_visits} is automatically added to the \code{fixed}
#'   formula.
#' @param method Character string specifying the estimation method passed
#'   to \code{\link[nlme]{lme}}. One of \code{"REML"} (restricted maximum
#'   likelihood; default) or \code{"ML"} (maximum likelihood).
#' @param optimizer Character string specifying the optimiser for the
#'   \code{\link[nlme]{lmeControl}} object. One of \code{"nlminb"}
#'   (default) or \code{"optim"}.
#'
#' @return A numeric matrix of fixed-effect estimates, identical to the
#'   \code{tTable} component of \code{summary(\link[nlme]{lme}(...))}. The
#'   matrix has one row per fixed effect and columns \code{Value},
#'   \code{Std.Error}, \code{DF}, \code{t-value}, and \code{p-value}.
#'
#' @details
#' This function is a convenience wrapper that:
#' \enumerate{
#'   \item Optionally computes a within-subject cumulative visit counter
#'     (\code{n_visits}) either before or after filtering by the
#'     observation indicator.
#'   \item Filters out rows where \code{obs_indicator == 0} (if
#'     specified).
#'   \item Removes incomplete cases.
#'   \item Fits a linear mixed model with \code{\link[nlme]{lme}} using
#'     the specified fixed and random formulas, estimation method, and
#'     optimiser.
#' }
#'
#' The OA-LMM and VA-LMM adjustments were proposed as simple corrections
#' for informative visit times. OA-LMM conditions on the number of
#' visits that yielded observed outcomes, while VA-LMM conditions on the
#' total number of visits. Neither fully accounts for the visiting
#' process but can reduce bias compared to a naive LMM.
#'
#' @references
#' Benjamin A. Goldstein, Nrupen A. Bhavsar, Matthew Phelan, and Michael J. Pencina (2016).
#' Controlling for informed presence bias due to the number of health encounters in an electronic health record.
#' \emph{American Journal of Epidemiology}, 184(11), 847--855. \doi{10.1093/aje/kww112}.
#'
#' Nan M. Laird and James H. Ware (1982). Random-effects models for longitudinal data.
#' \emph{Biometrics}, 38(4), 963--974.
#'
#' @examples
#' \donttest{
#' data(sim_ehr_data)
#' fit <- LMM(sim_ehr_data, fixed = Y ~ Z + X + time,
#'            random = ~ 1 | id, obs_indicator = "R")
#' print(fit)
#' summary(fit)
#' }
#'
#' @importFrom stats ave na.omit update
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

  result <- summary(fit)$tTable
  return(structure(list(coefficients = result), class = "LMM_result"))
}

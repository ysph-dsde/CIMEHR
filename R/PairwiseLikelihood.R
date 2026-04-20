#' Pairwise Composite Likelihood (PCL) Generalized Linear Model (GLM) Estimation
#'
#' @description
#' Fits a Generalized Linear Model (GLM) using Pairwise Composite
#' Likelihood (PCL). For each pair of subjects, the contribution to the
#' composite log-likelihood is based on the difference \eqn{Y_i - Y_j}
#' (Gaussian family) or a logistic-style softplus function (binomial
#' family). This approach yields consistent slope estimates without
#' requiring a correctly specified within-cluster correlation structure,
#' and is robust to informative cluster sizes.
#'
#' The intercept is not identifiable in pairwise differences and is
#' automatically removed if present in the formula.
#'
#' @param data A \code{data.frame} containing longitudinal or clustered
#'   data. Rows with missing values in the outcome, covariates, or
#'   identifier columns are removed.
#' @param formula A \code{\link[stats]{formula}} specifying the model
#'   (e.g., \code{Y ~ Z + X + time}). Must be supplied by the user;
#'   there is no default. If an intercept is present it is removed
#'   automatically with a message (when \code{verbose = TRUE}).
#' @param id_col Character string giving the name of the subject/cluster
#'   identifier column. Default is \code{"id"}.
#' @param time_col Character string giving the name of the time variable
#'   column. Default is \code{"time"}.
#' @param y_col Character string giving the name of the outcome column.
#'   Default is \code{"Y"}.
#' @param family Character string specifying the response distribution.
#'   One of:
#'   \describe{
#'     \item{\code{"gaussian"}}{(Default.) Pairwise residual sum of
#'       squares: \eqn{\sum 0.5\,(Y_i - Y_j - \mu_i + \mu_j)^2}.}
#'     \item{\code{"binomial"}}{Logistic-style softplus:
#'       \eqn{\sum -\log(1 + \exp((Y_i - Y_j)(\mu_j - \mu_i)))}.
#'       Requires \eqn{Y \in \{0, 1\}}.}
#'   }
#' @param start Optional numeric vector of starting values for the
#'   regression coefficients. Default is \code{NULL}, which uses Ordinary
#'   Least Squares (Gaussian) or GLM (binomial) estimates.
#' @param pair_sample Optional positive integer; if the total number of
#'   subject pairs exceeds this value, a random subsample of pairs is
#'   used. Default is \code{NULL} (use all pairs).
#' @param pair_seed Integer seed for reproducible pair subsampling.
#'   Default is \code{1}.
#' @param hessian Logical; if \code{TRUE}, the numeric Hessian is
#'   computed at the optimum via \code{\link[numDeriv]{hessian}} and used
#'   to obtain a naive variance--covariance matrix. Default is
#'   \code{TRUE}.
#' @param verbose Logical; if \code{TRUE}, progress messages are printed.
#'   Default is \code{FALSE}.
#'
#' @return An object of class \code{"PairwiseLikelihood"}, which is a
#'   named list with the following components:
#'
#' \describe{
#'   \item{\code{call}}{The matched call.}
#'   \item{\code{family}}{Character string echoing the response family.}
#'   \item{\code{formula}}{The formula used (intercept removed if
#'     applicable).}
#'   \item{\code{coef}}{Named numeric vector of estimated regression
#'     coefficients.}
#'   \item{\code{loglik}}{Numeric scalar; composite log-likelihood at the
#'     optimum (Gaussian: negated pairwise residual sum of squares;
#'     binomial: pairwise softplus sum).}
#'   \item{\code{converged}}{Logical; \code{TRUE} if
#'     \code{\link[stats]{optim}} converged.}
#'   \item{\code{hessian}}{(When \code{hessian = TRUE}.) Numeric matrix;
#'     Hessian of the objective at the optimum.}
#'   \item{\code{vcov_naive}}{(When \code{hessian = TRUE}.) Numeric
#'     matrix; inverse Hessian, a naive (model-based)
#'     variance--covariance estimator for \code{coef}.}
#'   \item{\code{n}}{Integer; number of observations used.}
#'   \item{\code{p}}{Integer; number of regression coefficients.}
#'   \item{\code{G}}{Integer; number of unique subjects/clusters.}
#'   \item{\code{k_pairs}}{Integer; number of subject pairs used.}
#'   \item{\code{pair_sample}}{Integer or \code{NULL}; echoing the
#'     subsampling parameter.}
#' }
#'
#' @details
#' The pairwise composite likelihood objective aggregates contributions
#' from all (or a subsample of) between-subject pairs. For each pair of
#' subjects \eqn{(i, j)} and all cross-products of their observations,
#' the contribution is:
#'
#' \strong{Gaussian:}
#' \deqn{l_{ij}(\beta) = -\frac{1}{2} \sum_{s \in i} \sum_{t \in j}
#'   \bigl[(Y_s - Y_t) - (\mu_s - \mu_t)\bigr]^2,}
#' where \eqn{\mu_s = X_s^\top \beta}. Minimising the total gives the
#' PCL estimator.
#'
#' \strong{Binomial:}
#' \deqn{l_{ij}(\beta) = \sum_{s \in i} \sum_{t \in j}
#'   -\log\bigl[1 + \exp\bigl((Y_s - Y_t)(\mu_t - \mu_s)\bigr)\bigr].}
#'
#' The inner loops are implemented in C++ via \pkg{Rcpp} for efficiency.
#' Optimization is performed with \code{\link[stats]{optim}} using BFGS.
#'
#' S3 methods \code{print} and \code{summary} are available for objects
#' of class \code{"PairwiseLikelihood"}.
#'
#' @note
#' This function requires the \pkg{Rcpp} and \pkg{numDeriv} packages. The
#' inline C++ functions \code{pw_total_gaussian_rcpp} and
#' \code{pw_total_softplus_rcpp} are compiled on first use.
#'
#' @references
#' Weining Shen, Suyu Liu, Yong Chen, and Jing Ning (2019). Regression analysis of
#' longitudinal data with outcome-dependent sampling and informative censoring.
#' \emph{Scandinavian Journal of Statistics}, 46(3), 831--847. \doi{doi/10.1111/sjos.12373}.
#'
#' Yong Chen, Jing Ning, and Chunyan Cai (2015). Regression analysis of longitudinal data
#' with irregular and informative observation times. \emph{Biostatistics},
#' 16(4), 727--739. \doi{10.1093/biostatistics/kxv008}.
#'
#' @examples
#' dat <- sim_ehr_data[sim_ehr_data$id <= 100 &
#'                     sim_ehr_data$R == 1 & !is.na(sim_ehr_data$R), ]
#' fit <- PairwiseLikelihood(dat, formula = Y ~ Z + X + time,
#'                           pair_sample = 200)
#' print(fit)
#'
#' @importFrom stats coef complete.cases glm.fit lm.fit model.frame
#'   model.matrix optim terms update
#' @export
PairwiseLikelihood <- function(
    data,
    formula,
    id_col      = "id",
    time_col    = "time",
    y_col       = "Y",
    family      = c("gaussian", "binomial"),
    start       = NULL,
    pair_sample = NULL,
    pair_seed   = 1L,
    hessian     = TRUE,
    verbose     = FALSE
){
  family <- match.arg(family)

  # 1) keep complete rows for variables used by formula + id/time/y
  df <- data
  needed <- all.vars(stats::terms(formula))
  stopifnot(all(c(id_col, time_col, y_col, needed) %in% names(df)))
  keep <- stats::complete.cases(df[, c(id_col, time_col, y_col, needed), drop = FALSE])
  df <- df[keep, , drop = FALSE]

  # 2) model matrix (drop intercept if present)
  has_intercept <- attr(stats::terms(formula), "intercept") == 1
  if (has_intercept) {
    if (verbose) message("Dropping intercept (not identifiable in pairwise objective).")
    formula <- stats::update(formula, . ~ . - 1)
  }
  # rebuild model frame/matrix AFTER possible update
  mf <- stats::model.frame(formula, data = df, na.action = stats::na.fail)
  X  <- stats::model.matrix(attr(mf, "terms"), data = mf)

  y  <- df[[y_col]]
  id <- df[[id_col]]

  # binomial: coerce y to {0,1}
  if (family == "binomial") {
    if (is.factor(y)) y <- as.numeric(y) - 1L
    y <- as.numeric(y)
    if (!all(y %in% c(0,1))) stop("For binomial family, y must be 0/1.")
  }

  # 3) per-subject indices
  groups_idx <- split(seq_len(nrow(df)), id)
  groups_idx <- Filter(function(ix) length(ix) > 0, groups_idx)
  G <- length(groups_idx)
  if (G < 2L) stop("Need at least two non-empty groups.")

  # Concatenate y and build 0-based offsets for C++
  y_concat <- numeric(0)
  offsets  <- integer(G + 1L)
  offsets[1] <- 0L
  for (g in seq_len(G)){
    idx <- groups_idx[[g]]
    y_concat <- c(y_concat, y[idx])
    offsets[g + 1L] <- offsets[g] + length(idx)
  }

  # Design per subject (for mu = X beta)
  Xi_list <- lapply(groups_idx, function(ix) X[ix, , drop = FALSE])

  # 4) subject-pairs
  subj_pairs <- utils::combn(seq_len(G), 2L)
  if (!is.null(pair_sample) && ncol(subj_pairs) > pair_sample) {
    set.seed(pair_seed)
    subj_pairs <- subj_pairs[, sample.int(ncol(subj_pairs), pair_sample), drop = FALSE]
  }
  g1 <- subj_pairs[1, ]
  g2 <- subj_pairs[2, ]

  # 5) start values
  if (is.null(start)) {
    b0 <- tryCatch({
      if (family == "gaussian") coef(stats::lm.fit(x = X, y = y))
      else coef(stats::glm.fit(x = X, y = y, family = stats::binomial()))
    }, error = function(e) rep(0, ncol(X)))
    b0[is.na(b0)] <- 0
  } else b0 <- start
  stopifnot(length(b0) == ncol(X))

  # helper: concat mu in subject order
  make_mu_concat <- function(beta){
    mus <- lapply(Xi_list, function(M) as.vector(M %*% beta))
    unlist(mus, use.names = FALSE)
  }

  # 6) objective
  if (family == "gaussian") {
    # Gaussian: minimize pairwise squared residuals
    obj <- function(beta){
      mu_concat <- make_mu_concat(beta)
      pw_total_gaussian_rcpp(y = y_concat,
                             mu = mu_concat,
                             offsets = offsets,
                             g1 = g1, g2 = g2)
    }
  } else { # binomial
    # Binomial: maximize pairwise softplus loglik (so minimize negative)
    obj <- function(beta){
      mu_concat <- make_mu_concat(beta)
      -pw_total_softplus_rcpp(y = y_concat,
                              mu = mu_concat,
                              offsets = offsets,
                              g1 = g1, g2 = g2)
    }
  }

  if (verbose) message(sprintf("[PairwiseLikelihood] family=%s, G=%d subjects, N=%d rows, pairs=%d",
                               family, G, nrow(df), ncol(subj_pairs)))

  # 7) optimize
  opt <- stats::optim(par = b0, fn = obj, method = "BFGS",
                      control = list(maxit = 200, reltol = 1e-8))
  bhat   <- opt$par
  names(bhat) <- colnames(X)
  loglik <- -opt$value
  conv   <- (opt$convergence == 0L)

  # 8) (optional) numeric Hessian and naive vcov
  H <- NULL; V_naive <- NULL
  if (hessian) {
    if (!requireNamespace("numDeriv", quietly = TRUE)) {
      stop("Package 'numDeriv' is required for Hessian computation. Install it with: install.packages('numDeriv')",
           call. = FALSE)
    }
    H <- tryCatch(numDeriv::hessian(obj, bhat), error = function(e) NULL)
    if (!is.null(H)) {
      rownames(H) <- colnames(H) <- colnames(X)
      V_naive <- tryCatch(solve(H), error = function(e) NULL)
      if (!is.null(V_naive)) rownames(V_naive) <- colnames(V_naive) <- colnames(X)
    }
  }

  structure(list(
    call        = match.call(),
    family      = family,
    formula     = formula,
    coef        = bhat,
    loglik      = loglik,
    converged   = conv,
    hessian     = H,
    vcov_naive  = V_naive,
    n           = nrow(df),
    p           = ncol(X),
    G           = G,
    k_pairs     = ncol(subj_pairs),
    pair_sample = pair_sample
  ), class = "PairwiseLikelihood")
}

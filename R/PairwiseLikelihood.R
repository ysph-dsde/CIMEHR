if (!requireNamespace("Rcpp", quietly = TRUE)) install.packages("Rcpp")
library(Rcpp)

Rcpp::cppFunction(code = '
#include <Rcpp.h>
using namespace Rcpp;

// numerically stable softplus: log(1 + exp(x))
inline double softplus(double x){
  double ax = std::fabs(x);
  return (x > 0.0 ? x : 0.0) + std::log1p(std::exp(-ax));
}

// Logistic-style pairwise loglik:
// sum_{pairs} -softplus( (yi - yj) * (muj - mui) )
// [[Rcpp::export]]
double pw_total_softplus_rcpp(const NumericVector& y,
                              const NumericVector& mu,
                              const IntegerVector& offsets,
                              const IntegerVector& g1,
                              const IntegerVector& g2){
  const int M = g1.size();
  double total = 0.0;

  for (int k = 0; k < M; ++k){
    const int si = g1[k] - 1;   // 0-based subject index
    const int sj = g2[k] - 1;

    const int ai0 = offsets[si];
    const int ai1 = offsets[si + 1];
    const int aj0 = offsets[sj];
    const int aj1 = offsets[sj + 1];

    for (int ii = ai0; ii < ai1; ++ii){
      const double yi  = y[ii];
      const double mui = mu[ii];
      for (int jj = aj0; jj < aj1; ++jj){
        const double z = (yi - y[jj]) * (mu[jj] - mui);
        total += -softplus(z);
      }
    }
  }
  return total;
}
')

Rcpp::cppFunction(code = '
#include <Rcpp.h>
using namespace Rcpp;

// Gaussian pairwise "likelihood":
// sum_{pairs} 0.5 * ( (Yi - Yj) - (mui - muj) )^2
// => minimize this total: pairwise residual sum of squares
// [[Rcpp::export]]
double pw_total_gaussian_rcpp(const NumericVector& y,
                              const NumericVector& mu,
                              const IntegerVector& offsets,
                              const IntegerVector& g1,
                              const IntegerVector& g2){
  const int M = g1.size();
  double total = 0.0;

  for (int k = 0; k < M; ++k){
    const int si = g1[k] - 1;   // 0-based
    const int sj = g2[k] - 1;

    const int ai0 = offsets[si];
    const int ai1 = offsets[si + 1];
    const int aj0 = offsets[sj];
    const int aj1 = offsets[sj + 1];

    for (int ii = ai0; ii < ai1; ++ii){
      const double yi  = y[ii];
      const double mui = mu[ii];
      for (int jj = aj0; jj < aj1; ++jj){
        const double diffY  = yi        - y[jj];
        const double diffMu = mui       - mu[jj];
        const double resid  = diffY     - diffMu;
        total += 0.5 * resid * resid;   // 0.5 * (resid^2)
      }
    }
  }
  return total;
}
')

#' Pairwise Composite Likelihood GLM Estimation
#'
#' Fits a generalized linear model using pairwise composite likelihood.
#'
#' @param data Data frame containing longitudinal/clustered data.
#' @param formula Model formula (e.g., Y ~ Z + X + time). Intercept is automatically
#'   removed as it is not identifiable in pairwise differences.
#' @param id_col Subject/cluster ID column name. Default is "id".
#' @param time_col Time variable column name. Default is "time".
#' @param y_col Outcome variable column name. Default is "Y".
#' @param family "gaussian" (continuous outcome) or "binomial" (binary outcome). Default is "gaussian".
#' @param start Optional starting values for coefficients. Default is NULL (uses lm/glm estimates).
#' @param pair_sample Optional integer to sample subject pairs for large datasets. Default is NULL (use all).
#' @param pair_seed Seed for reproducible pair sampling. Default is 1.
#' @param hessian If TRUE, compute Hessian and naive vcov. Default is TRUE.
#' @param verbose If TRUE, print progress messages. Default is FALSE.
#'
#' @return A list containing: coef, vcov_naive, loglik, converged, n, G, k_pairs, etc.
#'
#' @export
PairwiseLikelihood <- function(
    data,
    formula     = Y ~ Z + X + time,
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
    if (!requireNamespace("numDeriv", quietly = TRUE)) install.packages("numDeriv")
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

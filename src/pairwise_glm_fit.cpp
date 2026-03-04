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

// Gaussian pairwise "likelihood":
// sum_{pairs} 0.5 * ( (Yi - Yj) - (mui - muj) )^2
// [[Rcpp::export]]
double pw_total_gaussian_rcpp(const NumericVector& y,
                              const NumericVector& mu,
                              const IntegerVector& offsets,
                              const IntegerVector& g1,
                              const IntegerVector& g2){
  const int M = g1.size();
  double total = 0.0;

  for (int k = 0; k < M; ++k){
    const int si = g1[k] - 1;
    const int sj = g2[k] - 1;

    const int ai0 = offsets[si];
    const int ai1 = offsets[si + 1];
    const int aj0 = offsets[sj];
    const int aj1 = offsets[sj + 1];

    for (int ii = ai0; ii < ai1; ++ii){
      const double yi  = y[ii];
      const double mui = mu[ii];
      for (int jj = aj0; jj < aj1; ++jj){
        const double diffY  = yi  - y[jj];
        const double diffMu = mui - mu[jj];
        const double resid  = diffY - diffMu;
        total += 0.5 * resid * resid;
      }
    }
  }
  return total;
}

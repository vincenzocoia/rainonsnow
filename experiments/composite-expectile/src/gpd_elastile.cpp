// alpha-elastile of a three-parameter GPD.
//
// Root of  G(t) = (2 alpha / s) [ (1-2p) phi(t) + (1-p)(t - m) ]
//                + (1 - alpha) [ F(t) - p ],
// where the first bracket vanishes at the expectile and the second at the
// quantile, both increasing in t. So G is increasing, the root is unique, and it
// is bracketed by the quantile and the expectile -- which are supplied.
// Everything is closed form for the GPD.

#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector gpd_elastile_cpp(NumericVector p, double mu, double sigma, double xi,
                               double c_coef, double d_coef,
                               NumericVector blo, NumericVector bhi,
                               double tol = 1e-10, int maxit = 60) {
  int n = p.size();
  NumericVector out(n);
  if (xi >= 1.0 || sigma <= 0.0) { for (int i = 0; i < n; ++i) out[i] = NA_REAL; return out; }
  double m = mu + sigma / (1.0 - xi);
  bool z0 = std::fabs(xi) < 1e-12;

  auto surv = [&](double x) -> double {
    if (x < mu) return 1.0;
    if (z0) return std::exp(-(x - mu) / sigma);
    double u = 1.0 + xi * (x - mu) / sigma;
    return (u > 0.0) ? std::pow(u, -1.0 / xi) : 0.0;
  };
  auto dens = [&](double x) -> double {
    if (x < mu) return 0.0;
    if (z0) return std::exp(-(x - mu) / sigma) / sigma;
    double u = 1.0 + xi * (x - mu) / sigma;
    return (u > 0.0) ? std::pow(u, -1.0 / xi - 1.0) / sigma : 0.0;
  };
  auto phi = [&](double x) -> double {
    if (x < mu) return m - x;
    if (z0) return sigma * std::exp(-(x - mu) / sigma);
    double u = 1.0 + xi * (x - mu) / sigma;
    if (u <= 0.0) return 0.0;
    return (sigma + xi * (x - mu)) * std::pow(u, -1.0 / xi) / (1.0 - xi);
  };

  for (int i = 0; i < n; ++i) {
    double pp = p[i];
    double lo = blo[i], hi = bhi[i];
    if (ISNA(pp) || !R_finite(lo) || !R_finite(hi)) {
      out[i] = R_finite(lo) ? lo : hi; continue;
    }
    if (lo > hi) { double t = lo; lo = hi; hi = t; }
    if (hi - lo <= tol * std::max(1.0, std::fabs(hi))) { out[i] = 0.5 * (lo + hi); continue; }
    auto G = [&](double x) {
      return c_coef * ((1.0 - 2.0 * pp) * phi(x) + (1.0 - pp) * (x - m)) +
             d_coef * ((1.0 - surv(x)) - pp);
    };
    auto dG = [&](double x) {
      return c_coef * ((2.0 * pp - 1.0) * surv(x) + (1.0 - pp)) + d_coef * dens(x);
    };
    double x = 0.5 * (lo + hi);
    for (int it = 0; it < maxit; ++it) {
      double g = G(x);
      if (g < 0.0) lo = x; else hi = x;
      double sl = dG(x);
      double nx = (sl > 0.0) ? x - g / sl : 0.5 * (lo + hi);
      if (!R_finite(nx) || nx <= lo || nx >= hi) nx = 0.5 * (lo + hi);
      double d = std::fabs(nx - x); x = nx;
      if (d < tol * std::max(1.0, std::fabs(x))) break;
    }
    out[i] = x;
  }
  return out;
}

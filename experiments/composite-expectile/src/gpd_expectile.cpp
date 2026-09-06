// GPD expectile, in closed form throughout.
//
// For the three-parameter GPD with S(x) = (1 + xi (x-mu)/sigma)^(-1/xi),
//   phi(x) = E[(X - x)^+] = (sigma + xi (x - mu)) S(x) / (1 - xi)   for x >= mu,
//                         = m - x                                   for x <  mu,
// so the expectile identification equation k phi(x) + m - x = 0 needs no
// special functions at all. `anchor` replaces the model's mean m in the
// "+ m - x" term for the mean-anchored variant, leaving phi (the tail-local
// part) as the model's.

#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector gpd_expectile_cpp(NumericVector p, double mu, double sigma, double xi,
                                double anchor, double tol = 1e-10, int maxit = 60) {
  int n = p.size();
  NumericVector out(n);
  if (xi >= 1.0 || sigma <= 0.0) { for (int i = 0; i < n; ++i) out[i] = NA_REAL; return out; }
  double m = mu + sigma / (1.0 - xi);
  bool zero_xi = std::fabs(xi) < 1e-12;
  double upper = (xi < 0) ? mu - sigma / xi : R_PosInf;

  auto surv = [&](double x) -> double {
    if (x < mu) return 1.0;
    if (zero_xi) return std::exp(-(x - mu) / sigma);
    double u = 1.0 + xi * (x - mu) / sigma;
    return (u > 0.0) ? std::pow(u, -1.0 / xi) : 0.0;
  };
  auto phi = [&](double x) -> double {
    if (x < mu) return m - x;
    if (zero_xi) return sigma * std::exp(-(x - mu) / sigma);
    double u = 1.0 + xi * (x - mu) / sigma;
    if (u <= 0.0) return 0.0;
    return (sigma + xi * (x - mu)) * std::pow(u, -1.0 / xi) / (1.0 - xi);
  };
  auto quant = [&](double pp) -> double {
    if (zero_xi) return mu - sigma * std::log1p(-pp);
    return mu + sigma * (std::pow(1.0 - pp, -xi) - 1.0) / xi;
  };

  double step = std::max(1.0, std::max(std::fabs(anchor), sigma));
  for (int i = 0; i < n; ++i) {
    double pp = p[i];
    if (ISNA(pp) || ISNAN(pp)) { out[i] = NA_REAL; continue; }
    if (pp <= 0.0) { out[i] = mu; continue; }
    if (pp >= 1.0) { out[i] = upper; continue; }
    double k = (2.0 * pp - 1.0) / (1.0 - pp);
    auto g = [&](double x) { return k * phi(x) + anchor - x; };

    double lo = anchor, hi = std::max(anchor + step, quant(pp) + step);
    if (g(lo) < 0.0) {                      // root below the anchor
      hi = anchor; lo = anchor - step;
      for (int it = 0; it < 300 && g(lo) < 0.0; ++it) {
        lo = anchor - (anchor - lo) * 4.0; if (lo < -1e300) break; }
    } else {
      for (int it = 0; it < 300 && g(hi) > 0.0; ++it) {
        hi = anchor + (hi - anchor) * 4.0; if (hi > 1e300) break; }
    }
    double x = quant(pp);
    if (!(x > lo && x < hi)) x = 0.5 * (lo + hi);
    for (int it = 0; it < maxit; ++it) {
      double gx = g(x);
      if (gx > 0.0) lo = x; else hi = x;
      double dg = -k * surv(x) - 1.0;
      double nx = x - gx / dg;
      if (!R_finite(nx) || nx <= lo || nx >= hi) nx = 0.5 * (lo + hi);
      double d = std::fabs(nx - x);
      x = nx;
      if (d < tol * std::max(1.0, std::fabs(x))) break;
    }
    out[i] = x;
  }
  return out;
}

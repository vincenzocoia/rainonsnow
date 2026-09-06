// Mean-anchored GEV expectile.
//
// The p-expectile solves  k phi(x) + m - x = 0  with m the distribution's mean.
// In the composite criterion that mean is the MODEL's, and under body
// contamination the model cannot reproduce it -- which is what leaves a
// residual bias in the fitted tail at every level (the "mean anchoring"
// problem). Here m is replaced by an externally supplied anchor, so that the
// tail-local part (phi) comes from the model, where the GEV is right, and the
// body-dominated part (the mean) comes from the data, where it is not.
//
// In standard coordinates x = mu + sigma tau the equation is
//   k phi_0(tau) + mtilde - tau = 0,   mtilde = (m_anchor - mu) / sigma,
// so the solver is the ordinary one with the mean replaced by mtilde. It is
// still strictly decreasing in tau, so the same safeguarded Newton applies.

#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

static inline double zc(double x, double xi) {
  double t = 1.0 + xi * x;
  if (t <= 0.0) return (xi > 0.0) ? R_PosInf : 0.0;
  return std::pow(t, -1.0 / xi);
}
static inline double phic(double x, double xi, double m0, double lo, double hi, double gam) {
  if (x <= lo) return m0 - x;
  if (x >= hi) return 0.0;
  double z = zc(x, xi);
  if (!R_finite(z)) return m0 - x;
  return (R::pgamma(z, 1.0 - xi, 1.0, 1, 0) * gam - (-std::expm1(-z)) * std::pow(z, -xi)) / xi;
}

// [[Rcpp::export]]
NumericVector gev_expectile_anchored_cpp(NumericVector p, double xi, double mtilde,
                                         double tol = 1e-11, int maxit = 100) {
  int n = p.size();
  NumericVector out(n);
  if (std::fabs(xi) < 1e-7) xi = (xi >= 0 ? 1e-7 : -1e-7);
  if (xi >= 1.0) { for (int i = 0; i < n; ++i) out[i] = NA_REAL; return out; }
  double gam = ::tgamma(1.0 - xi);
  double m0 = (gam - 1.0) / xi;                 // the model's own standard mean
  double lo_end = (xi > 0) ? -1.0 / xi : R_NegInf;
  double hi_end = (xi < 0) ? -1.0 / xi : R_PosInf;

  for (int i = 0; i < n; ++i) {
    double pp = p[i];
    if (ISNA(pp) || ISNAN(pp)) { out[i] = NA_REAL; continue; }
    if (pp <= 0.0) { out[i] = lo_end; continue; }
    if (pp >= 1.0) { out[i] = hi_end; continue; }
    double k = (2.0 * pp - 1.0) / (1.0 - pp);
    // g(x) = k phi(x) + mtilde - x, strictly decreasing.
    auto g = [&](double x) { return k * phic(x, xi, m0, lo_end, hi_end, gam) + mtilde - x; };
    double lo, hi;
    // g(mtilde) = k phi(mtilde) >= 0, so the root is at or above the anchor.
    lo = mtilde; hi = mtilde + 1.0;
    for (int it = 0; it < 300; ++it) { if (g(hi) <= 0.0) break;
      hi = mtilde + (hi - mtilde) * 4.0; if (hi > 1e300) break; }
    if (g(lo) < 0.0) {                          // p below 1/2 can put it lower
      hi = mtilde; lo = mtilde - 1.0;
      for (int it = 0; it < 300; ++it) { if (g(lo) >= 0.0) break;
        lo = mtilde - (mtilde - lo) * 4.0; if (lo < -1e300) break; }
      if (R_finite(lo_end) && lo < lo_end) lo = lo_end;
    }
    double x = 0.5 * (lo + hi);
    if (!R_finite(x)) x = mtilde;
    for (int it = 0; it < maxit; ++it) {
      double gx = g(x);
      if (gx > 0.0) lo = x; else hi = x;
      double dg = -k * (-std::expm1(-zc(x, xi))) - 1.0;
      double nx = x - gx / dg;
      if (!R_finite(nx) || nx <= lo || nx >= hi) nx = 0.5 * (lo + hi);
      double step = std::fabs(nx - x);
      x = nx;
      if (step < tol * std::max(1.0, std::fabs(x))) break;
    }
    out[i] = x;
  }
  return out;
}

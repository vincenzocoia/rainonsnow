// The alpha-elastile of a GEV: the functional elicited by a convex combination
// of the expectile (L2) and quantile (L1) check functions at the same level.
//
// At level p the combined loss of a prediction t is
//
//   rho^alpha_p(y, t) = (alpha/s) |p - I(y<t)| (y - t)^2
//                     + (1 - alpha) (p - I(y<t)) (y - t),
//
// where s carries units of y so the two terms are commensurable. Differentiating
// the expected loss and writing phi(t) = E[(Y-t)^+], m = E[Y]:
//
//   G(t) = (2 alpha / s) [ (1 - 2p) phi(t) + (1 - p)(t - m) ]
//        + (1 - alpha)   [ F(t) - p ]                          = 0.
//
// The first bracket vanishes at the p-expectile, the second at the p-quantile,
// and each is increasing in t. So G is strictly increasing -- the root is unique
// -- and it is bracketed by the quantile and the expectile themselves, which is
// the bracket used below. alpha = 0 gives the quantile, alpha = 1 the expectile.
//
// Location-scale: with t = mu + sigma tau the equation becomes
//   c A_0(tau) + d B_0(tau) = 0,  c = 2 alpha sigma / s,  d = 1 - alpha,
// in the standard GEV(0, 1, xi), which is what this solves.

#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

static inline double z_of(double x, double xi) {
  double t = 1.0 + xi * x;
  if (t <= 0.0) return (xi > 0.0) ? R_PosInf : 0.0;
  return std::pow(t, -1.0 / xi);
}

static inline double phi_std(double x, double xi, double m, double lo_end,
                             double hi_end, double gam) {
  if (x <= lo_end) return m - x;
  if (x >= hi_end) return 0.0;
  double z = z_of(x, xi);
  if (!R_finite(z)) return m - x;
  return (R::pgamma(z, 1.0 - xi, 1.0, 1, 0) * gam -
          (-std::expm1(-z)) * std::pow(z, -xi)) / xi;
}

// [[Rcpp::export]]
NumericVector gev_elastile_std_cpp(NumericVector p, double xi, double c, double d,
                                   NumericVector bracket_lo, NumericVector bracket_hi,
                                   double tol = 1e-11, int maxit = 100) {
  int n = p.size();
  NumericVector out(n);
  if (std::fabs(xi) < 1e-7) xi = (xi >= 0 ? 1e-7 : -1e-7);
  if (xi >= 1.0) { for (int i = 0; i < n; ++i) out[i] = NA_REAL; return out; }

  double gam = ::tgamma(1.0 - xi);
  double m = (gam - 1.0) / xi;
  double lo_end = (xi > 0) ? -1.0 / xi : R_NegInf;
  double hi_end = (xi < 0) ? -1.0 / xi : R_PosInf;

  for (int i = 0; i < n; ++i) {
    double pp = p[i];
    double lo = bracket_lo[i], hi = bracket_hi[i];
    if (ISNA(pp) || ISNAN(pp) || !R_finite(lo) || !R_finite(hi)) {
      // A non-finite bracket end means the corresponding functional is at an
      // infinite endpoint; fall back to whichever end is finite.
      out[i] = R_finite(lo) ? lo : hi;
      continue;
    }
    if (lo > hi) { double t = lo; lo = hi; hi = t; }
    if (hi - lo <= tol * std::max(1.0, std::fabs(hi))) { out[i] = 0.5 * (lo + hi); continue; }

    auto G = [&](double x) {
      double z = z_of(x, xi);
      double F = std::exp(-z);
      return c * ((1.0 - 2.0 * pp) * phi_std(x, xi, m, lo_end, hi_end, gam) +
                  (1.0 - pp) * (x - m)) + d * (F - pp);
    };
    auto dG = [&](double x) {
      double z = z_of(x, xi);
      double S = -std::expm1(-z);
      // density of the standard GEV: z^(xi+1) e^{-z}
      double f = (R_finite(z) && z > 0) ? std::pow(z, xi + 1.0) * std::exp(-z) : 0.0;
      return c * ((2.0 * pp - 1.0) * S + (1.0 - pp)) + d * f;
    };

    double x = 0.5 * (lo + hi);
    for (int it = 0; it < maxit; ++it) {
      double g = G(x);
      if (g < 0.0) lo = x; else hi = x;          // G is increasing
      double slope = dG(x);
      double nx = (slope > 0.0) ? x - g / slope : 0.5 * (lo + hi);
      if (!R_finite(nx) || nx <= lo || nx >= hi) nx = 0.5 * (lo + hi);
      double step = std::fabs(nx - x);
      x = nx;
      if (step < tol * std::max(1.0, std::fabs(x))) break;
    }
    out[i] = x;
  }
  return out;
}

// Fast GEV expectile function.
//
// Mirrors gev_expectile_std() in R/gev.R exactly: the p-expectile of the
// standard GEV(0, 1, xi) is the unique root in x of
//
//   g_p(x) = k phi(x) + m - x,     k = (2p - 1) / (1 - p),
//
// with phi(x) = E[(X - x)^+] in closed form via the lower incomplete gamma,
//
//   phi(x) = (1/xi) [ gammainc_lower(1 - xi, z) - (1 - e^{-z}) z^{-xi} ],
//   z = (1 + xi x)^(-1/xi).
//
// g_p is strictly decreasing with g_p'(x) = -k S(x) - 1, so a safeguarded
// Newton iteration (bisect whenever the Newton step leaves the bracket)
// converges from any bracket. Location and scale are applied on the R side.

#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

static inline double gev_z_c(double x, double xi) {
  double t = 1.0 + xi * x;
  if (t <= 0.0) return (xi > 0.0) ? R_PosInf : 0.0;
  return std::pow(t, -1.0 / xi);
}

static inline double phi_c(double x, double xi, double m,
                           double lo_end, double hi_end, double gam1mxi) {
  if (x <= lo_end) return m - x;    // X always exceeds x
  if (x >= hi_end) return 0.0;      // X never exceeds x
  double z = gev_z_c(x, xi);
  if (!R_finite(z)) return m - x;
  double gl = R::pgamma(z, 1.0 - xi, 1.0, 1, 0) * gam1mxi;
  return (gl - (-std::expm1(-z)) * std::pow(z, -xi)) / xi;
}

static inline double surv_c(double x, double xi) {
  return -std::expm1(-gev_z_c(x, xi));
}

// [[Rcpp::export]]
NumericVector gev_expectile_std_cpp(NumericVector p, double xi,
                                    double tol = 1e-11, int maxit = 100) {
  int n = p.size();
  NumericVector out(n);
  // The closed form for phi cancels badly as xi -> 0; nudging xi by 1e-7 costs
  // about nine significant digits' headroom and is far below any scale that
  // matters here.
  if (std::fabs(xi) < 1e-7) xi = (xi >= 0 ? 1e-7 : -1e-7);
  if (xi >= 1.0) { for (int i = 0; i < n; ++i) out[i] = NA_REAL; return out; }

  double gam1mxi = ::tgamma(1.0 - xi);
  double m = (gam1mxi - 1.0) / xi;
  double lo_end = (xi > 0) ? -1.0 / xi : R_NegInf;
  double hi_end = (xi < 0) ? -1.0 / xi : R_PosInf;

  for (int i = 0; i < n; ++i) {
    double pp = p[i];
    if (ISNA(pp) || ISNAN(pp)) { out[i] = NA_REAL; continue; }
    if (pp <= 0.0) { out[i] = lo_end; continue; }
    if (pp >= 1.0) { out[i] = hi_end; continue; }
    if (pp == 0.5) { out[i] = m; continue; }

    double k = (2.0 * pp - 1.0) / (1.0 - pp);
    double lo, hi;
    if (pp > 0.5) {
      lo = m; hi = m + 1.0;
      // Expand upward until g < 0. The root is above the mean for p > 1/2.
      for (int it = 0; it < 200; ++it) {
        double g = k * phi_c(hi, xi, m, lo_end, hi_end, gam1mxi) + m - hi;
        if (g <= 0.0) break;
        hi = m + (hi - m) * 4.0;
        if (hi > 1e300) break;
      }
    } else {
      hi = m; lo = m - 1.0;
      for (int it = 0; it < 200; ++it) {
        if (R_finite(lo_end) && lo <= lo_end) { lo = lo_end; break; }
        double g = k * phi_c(lo, xi, m, lo_end, hi_end, gam1mxi) + m - lo;
        if (g >= 0.0) break;
        lo = m - (m - lo) * 4.0;
        if (lo < -1e300) break;
      }
      if (R_finite(lo_end) && lo < lo_end) lo = lo_end;
    }

    double x = 0.5 * (lo + hi);
    if (!R_finite(x)) x = m;
    for (int it = 0; it < maxit; ++it) {
      double g = k * phi_c(x, xi, m, lo_end, hi_end, gam1mxi) + m - x;
      if (g > 0.0) lo = x; else hi = x;
      double dg = -k * surv_c(x, xi) - 1.0;
      double nx = x - g / dg;
      if (!R_finite(nx) || nx <= lo || nx >= hi) nx = 0.5 * (lo + hi);
      double step = std::fabs(nx - x);
      x = nx;
      if (step < tol * std::max(1.0, std::fabs(x))) break;
    }
    out[i] = x;
  }
  return out;
}

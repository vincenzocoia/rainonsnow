# ---------------------------------------------------------------------------
# GEV distribution: cdf, quantile, mean, partial moment, expectile.
#
# Parameterisation
#   F(x) = exp(-(1 + xi (x - mu) / sigma)^(-1/xi)),   1 + xi (x - mu)/sigma > 0,
# with the Gumbel limit at xi = 0.
#
# The expectile function is needed on a grid of levels inside an optimiser, so
# the routines here are vectorised over their first argument and the expectile
# solver is a vectorised safeguarded Newton iteration.
# ---------------------------------------------------------------------------

EULER <- 0.577215664901532861

gev_z <- function(x, mu, sigma, xi) {
  # z = (1 + xi (x - mu)/sigma)^(-1/xi), so that F(x) = exp(-z).
  # z decreases from +Inf at the lower endpoint to 0 at the upper endpoint.
  if (abs(xi) < 1e-12) return(exp(-(x - mu) / sigma))
  t <- 1 + xi * (x - mu) / sigma
  z <- rep(NA_real_, length(t))
  ok <- t > 0
  z[ok] <- t[ok]^(-1 / xi)
  z[!ok] <- if (xi > 0) Inf else 0
  z
}

pgev <- function(x, mu, sigma, xi) exp(-gev_z(x, mu, sigma, xi))
gev_survival <- function(x, mu, sigma, xi) -expm1(-gev_z(x, mu, sigma, xi))

qgev <- function(p, mu, sigma, xi) {
  if (abs(xi) < 1e-12) return(mu - sigma * log(-log(p)))
  mu + sigma * ((-log(p))^(-xi) - 1) / xi
}

gev_mean <- function(mu, sigma, xi) {
  if (xi >= 1) return(Inf)
  if (abs(xi) < 1e-12) return(mu + sigma * EULER)
  mu + sigma * (gamma(1 - xi) - 1) / xi
}

gev_endpoints <- function(mu, sigma, xi) {
  if (xi > 1e-12) c(mu - sigma / xi, Inf)
  else if (xi < -1e-12) c(-Inf, mu - sigma / xi)
  else c(-Inf, Inf)
}

# ---------------------------------------------------------------------------
# Partial moment  phi(x) = E[(X - x)^+] = \int_x^Inf S(t) dt.
#
# Substituting z = (1 + xi (t - mu)/sigma)^(-1/xi)  (dt = -sigma z^{-xi-1} dz)
# gives  phi(x) = sigma \int_0^{z_x} (1 - e^{-z}) z^{-xi-1} dz, and integrating
# by parts (u = 1 - e^{-z}, dv = z^{-xi-1} dz) gives the closed form
#
#   phi(x) = (sigma / xi) [ gammainc_lower(1 - xi, z_x) - (1 - e^{-z_x}) z_x^{-xi} ]
#
# for xi < 1, xi != 0. The boundary term vanishes at z = 0 because
# (1 - e^{-z}) z^{-xi} ~ z^{1-xi} -> 0. As x falls to the lower endpoint,
# z_x -> Inf and phi -> sigma Gamma(1 - xi) / xi = mean - lower endpoint.
#
# At xi = 0 (Gumbel) the same integral is
#   phi(x) = sigma [ log(z) + euler + E1(z) ]   with  E1 the exponential integral.
# ---------------------------------------------------------------------------
gev_partial_moment <- function(x, mu, sigma, xi) {
  stopifnot(xi < 1)
  ends <- gev_endpoints(mu, sigma, xi)
  m <- gev_mean(mu, sigma, xi)
  out <- numeric(length(x))
  below <- x <= ends[1]        # X is always above x, so E[(X-x)^+] = m - x
  above <- x >= ends[2]        # X is never above x
  mid <- !below & !above
  out[below] <- m - x[below]
  out[above] <- 0
  if (any(mid)) {
    z <- gev_z(x[mid], mu, sigma, xi)
    if (abs(xi) < 1e-8) {
      out[mid] <- sigma * (log(z) + EULER + expint_E1(z))
    } else {
      gl <- pgamma(z, shape = 1 - xi) * gamma(1 - xi)   # lower incomplete gamma
      out[mid] <- (sigma / xi) * (gl - (-expm1(-z)) * z^(-xi))
    }
  }
  out
}

# E1(a) = \int_a^Inf e^{-t}/t dt. Only used in the (measure-zero) Gumbel case.
expint_E1 <- function(a) {
  vapply(a, function(aa) {
    if (!is.finite(aa)) return(0)
    if (aa <= 0) return(Inf)
    integrate(function(t) exp(-t) / t, aa, Inf, rel.tol = 1e-10)$value
  }, numeric(1))
}

# ---------------------------------------------------------------------------
# Expectile function.
#
# The p-expectile is the unique root in x of the identification equation
#   g_p(x) = k phi(x) + m - x,      k = (2p - 1)/(1 - p),
# which is strictly decreasing with g_p'(x) = -k S(x) - 1.
#
# GEV is a location-scale family in (mu, sigma) for fixed xi, and expectiles are
# location-scale equivariant, so we solve once in the standard GEV(0, 1, xi) and
# rescale. That keeps the shape-dependent work in one place.
# ---------------------------------------------------------------------------
gev_expectile_std <- function(p, xi, tol = 1e-11, maxit = 60) {
  m <- gev_mean(0, 1, xi)
  if (!is.finite(m)) return(rep(NA_real_, length(p)))
  ends <- gev_endpoints(0, 1, xi)
  n <- length(p)
  out <- rep(NA_real_, n)
  ok <- !is.na(p) & p > 0 & p < 1
  out[!is.na(p) & p <= 0] <- ends[1]
  out[!is.na(p) & p >= 1] <- ends[2]
  if (!any(ok)) return(out)
  pp <- p[ok]
  k <- (2 * pp - 1) / (1 - pp)
  g  <- function(x) k * gev_partial_moment(x, 0, 1, xi) + m - x
  dg <- function(x) -k * gev_survival(x, 0, 1, xi) - 1

  # Bracket. g is decreasing in x, and x = m is the p = 1/2 solution, so the
  # root lies above m when p > 1/2 and below m when p < 1/2.
  lo <- ifelse(pp > 0.5, m, ends[1]); hi <- ifelse(pp > 0.5, Inf, m)
  # Finite upper bracket for the p > 1/2 side by geometric expansion.
  up <- rep(m + 1, sum(pp > 0.5))
  if (any(pp > 0.5)) {
    idx <- which(pp > 0.5)
    for (it in 1:200) {
      gv <- k[idx] * gev_partial_moment(up, 0, 1, xi) + m - up
      grow <- gv > 0
      if (!any(grow)) break
      up[grow] <- m + (up[grow] - m) * 4
      if (max(up) > 1e18) break
    }
    hi[idx] <- up
  }
  if (!is.finite(ends[1])) {
    idx <- which(pp <= 0.5)
    if (length(idx)) {
      dn <- rep(m - 1, length(idx))
      for (it in 1:200) {
        gv <- k[idx] * gev_partial_moment(dn, 0, 1, xi) + m - dn
        grow <- gv < 0
        if (!any(grow)) break
        dn[grow] <- m - (m - dn[grow]) * 4
        if (min(dn) < -1e18) break
      }
      lo[idx] <- dn
    }
  }
  x <- (lo + hi) / 2
  x[!is.finite(x)] <- m
  # Safeguarded Newton: take the Newton step when it stays inside the bracket,
  # otherwise bisect. Vectorised over levels.
  for (it in seq_len(maxit)) {
    gx <- g(x)
    pos <- gx > 0
    lo[pos] <- x[pos]; hi[!pos] <- x[!pos]
    nx <- x - gx / dg(x)
    bad <- !is.finite(nx) | nx <= lo | nx >= hi
    nx[bad] <- ((lo + hi) / 2)[bad]
    if (max(abs(nx - x) / pmax(1, abs(x))) < tol) { x <- nx; break }
    x <- nx
  }
  out[ok] <- x
  out
}

# The C++ solver (src/gev_expectile.cpp) is a line-for-line port of
# gev_expectile_std() above and is used when compiled; the R version is kept as
# the readable reference and is checked against it in scripts/98-validate.R.
gev_expectile <- function(p, mu, sigma, xi, ...) {
  std <- if (exists("gev_expectile_std_cpp", mode = "function"))
    gev_expectile_std_cpp(p, xi) else gev_expectile_std(p, xi, ...)
  mu + sigma * std
}

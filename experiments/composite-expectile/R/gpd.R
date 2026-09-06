# ---------------------------------------------------------------------------
# Generalised Pareto distribution, three-parameter (location mu, scale sigma,
# shape xi):
#
#   S(x) = (1 + xi (x - mu) / sigma)^(-1/xi),   x >= mu,
#
# with the exponential limit at xi = 0. Everything the composite estimators
# need is closed form here, which makes the GPD cheaper to work with than the
# GEV: no incomplete gamma, no numerical integration anywhere.
#
#   phi(x) = E[(X - x)^+]           = (sigma + xi (x - mu)) S(x) / (1 - xi)
#   psi(x) = E[(X - x)^2 1{X >= x}] = 2 sigma^2 u^2 S(x) / ((1-xi)(1-2 xi))
#
# with u = 1 + xi (x - mu)/sigma. Both are verified against numerical
# integration and Monte Carlo in scripts/98-validate.R.
# ---------------------------------------------------------------------------

gpd_u <- function(x, mu, sigma, xi) 1 + xi * (x - mu) / sigma

sgpd <- function(x, mu, sigma, xi) {          # survival
  out <- numeric(length(x))
  below <- x < mu
  out[below] <- 1
  ok <- !below
  if (any(ok)) {
    if (abs(xi) < 1e-12) {
      out[ok] <- exp(-(x[ok] - mu) / sigma)
    } else {
      u <- gpd_u(x[ok], mu, sigma, xi)
      out[ok] <- ifelse(u > 0, u^(-1 / xi), 0)
    }
  }
  out
}
pgpd <- function(x, mu, sigma, xi) 1 - sgpd(x, mu, sigma, xi)

qgpd <- function(p, mu, sigma, xi) {
  if (abs(xi) < 1e-12) return(mu - sigma * log(1 - p))
  mu + sigma * ((1 - p)^(-xi) - 1) / xi
}

dgpd <- function(x, mu, sigma, xi) {
  out <- numeric(length(x))
  ok <- x >= mu
  if (!any(ok)) return(out)
  if (abs(xi) < 1e-12) { out[ok] <- exp(-(x[ok] - mu) / sigma) / sigma; return(out) }
  u <- gpd_u(x[ok], mu, sigma, xi)
  out[ok] <- ifelse(u > 0, u^(-1 / xi - 1) / sigma, 0)
  out
}

gpd_mean <- function(mu, sigma, xi) if (xi >= 1) Inf else mu + sigma / (1 - xi)

gpd_endpoints <- function(mu, sigma, xi)
  c(mu, if (xi < 0) mu - sigma / xi else Inf)

gpd_partial_moment <- function(x, mu, sigma, xi) {
  stopifnot(xi < 1)
  m <- gpd_mean(mu, sigma, xi)
  out <- numeric(length(x))
  below <- x < mu
  out[below] <- m - x[below]
  ok <- !below
  if (any(ok)) {
    if (abs(xi) < 1e-12) {
      out[ok] <- sigma * exp(-(x[ok] - mu) / sigma)
    } else {
      out[ok] <- (sigma + xi * (x[ok] - mu)) * sgpd(x[ok], mu, sigma, xi) / (1 - xi)
    }
  }
  out
}

gpd_second_partial_moment <- function(x, mu, sigma, xi) {
  if (xi >= 0.5) return(rep(Inf, length(x)))
  m <- gpd_mean(mu, sigma, xi)
  v2 <- 2 * sigma^2 / ((1 - xi) * (1 - 2 * xi))     # E[(X - mu)^2]
  out <- numeric(length(x))
  below <- x < mu
  out[below] <- v2 - 2 * (x[below] - mu) * (m - mu) + (x[below] - mu)^2
  ok <- !below
  if (any(ok)) {
    u <- gpd_u(x[ok], mu, sigma, xi)
    out[ok] <- 2 * sigma^2 * u^2 * sgpd(x[ok], mu, sigma, xi) / ((1 - xi) * (1 - 2 * xi))
  }
  out
}

# Expectile: the unique root of k phi(x) + m - x = 0, k = (2p-1)/(1-p).
# Strictly decreasing, and with phi in closed form a vectorised safeguarded
# Newton needs only a handful of iterations. `m_anchor` replaces the model's
# mean for the mean-anchored variant.
# The R implementation below is the readable reference; src/gpd_expectile.cpp is
# a line-for-line port used when compiled, and the two are checked against each
# other in scripts/98-validate.R.
gpd_expectile <- function(p, mu, sigma, xi, m_anchor = NULL, tol = 1e-11, maxit = 60) {
  if (exists("gpd_expectile_cpp", mode = "function")) {
    anchor <- if (is.null(m_anchor)) gpd_mean(mu, sigma, xi) else m_anchor
    if (!is.finite(anchor)) return(rep(NA_real_, length(p)))
    return(gpd_expectile_cpp(p, mu, sigma, xi, anchor))
  }
  m <- gpd_mean(mu, sigma, xi)
  if (!is.finite(m)) return(rep(NA_real_, length(p)))
  anchor <- if (is.null(m_anchor)) m else m_anchor
  ends <- gpd_endpoints(mu, sigma, xi)
  out <- rep(NA_real_, length(p))
  ok <- !is.na(p) & p > 0 & p < 1
  out[!is.na(p) & p <= 0] <- ends[1]; out[!is.na(p) & p >= 1] <- ends[2]
  if (!any(ok)) return(out)
  pp <- p[ok]; k <- (2 * pp - 1) / (1 - pp)
  g <- function(x) k * gpd_partial_moment(x, mu, sigma, xi) + anchor - x
  step <- max(1, abs(anchor), sigma)
  # Bracket from the quantile rather than by blind expansion: the expectile and
  # the quantile at the same level are close, so this usually needs no
  # expansion at all, where a geometric search from the mean needs several
  # passes over the whole level grid.
  qq <- qgpd(pp, mu, sigma, xi)
  lo <- rep(anchor, length(pp))
  hi <- pmax(anchor + step, qq + step)
  for (it in 1:300) { bad <- g(hi) > 0; if (!any(bad)) break
    hi[bad] <- anchor + (hi[bad] - anchor) * 4; if (max(hi) > 1e300) break }
  low <- g(lo) < 0
  if (any(low)) {
    hi[low] <- anchor; lo[low] <- anchor - step
    for (it in 1:300) { bad <- low & (g(lo) < 0); if (!any(bad)) break
      lo[bad] <- anchor - (anchor - lo[bad]) * 4; if (min(lo) < -1e300) break }
    # No clipping to the support's lower end: below mu the partial moment is
    # m - x, so g stays linear and the root is still well defined. With an
    # anchor supplied from data the root can legitimately sit below mu.
  }
  x <- pmin(pmax(qq, lo + 1e-12), hi - 1e-12)      # start from the quantile
  x[!is.finite(x)] <- ((lo + hi) / 2)[!is.finite(x)]
  for (it in seq_len(maxit)) {
    gx <- g(x); pos <- gx > 0
    lo[pos] <- x[pos]; hi[!pos] <- x[!pos]
    dg <- -k * sgpd(x, mu, sigma, xi) - 1
    nx <- x - gx / dg
    bad <- !is.finite(nx) | nx <= lo | nx >= hi
    nx[bad] <- ((lo + hi) / 2)[bad]
    if (max(abs(nx - x) / pmax(1, abs(x))) < tol) { x <- nx; break }
    x <- nx
  }
  out[ok] <- x
  out
}

# alpha-elastile of the GPD; alpha = 0 is the quantile, 1 the expectile, and the
# root always lies between them. `s` makes the two loss terms commensurable.
gpd_elastile <- function(p, mu, sigma, xi, alpha, s) {
  if (alpha <= 0) return(qgpd(p, mu, sigma, xi))
  if (alpha >= 1) return(gpd_expectile(p, mu, sigma, xi))
  q <- qgpd(p, mu, sigma, xi); e <- gpd_expectile(p, mu, sigma, xi)
  gpd_elastile_cpp(p, mu, sigma, xi, 2 * alpha / s, 1 - alpha,
                   pmin(q, e), pmax(q, e))
}

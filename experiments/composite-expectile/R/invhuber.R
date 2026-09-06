# ---------------------------------------------------------------------------
# Inverted-Huber M-quantiles.
#
# Breckling & Chambers (1988) define the M-quantile of level p through
#
#     E[ |p - I(Y < t)| psi(Y - t) ] = 0
#
# for a chosen influence function psi. psi(u) = sign(u) gives the quantile,
# psi(u) = u the expectile, and Huber's psi -- linear near zero, bounded beyond
# a knot c -- gives the M-quantile proper. Huber exists to bound the influence
# of large residuals. In tail estimation that is backwards: the large
# observations carry nearly all the information about xi, and bounding them
# discards exactly what is wanted. So invert the knot.
#
# Two distinct inversions are possible, and they are not the same estimator.
#
# (A) SYMMETRIC, psi(u) = sign(u) max(|u|, c).
#     Bounded (constant) influence within c of the fitted level, linear beyond
#     it on both sides. Continuous at |u| = c. Endpoints: c -> 0 is the
#     expectile, c -> Inf the quantile -- the same two ends as Huber, reached in
#     the opposite order and by a different route in between.
#
# (B) ONE-SIDED, psi(u) = max(u, -c).
#     Linear above the fitted level, capped at -c below it. This is the one the
#     "bounded in the body, unbounded in the tail" argument actually asks for:
#     it is the downside that carries the contamination. Endpoints: c -> Inf is
#     the expectile; c -> 0 is degenerate (see below).
#
# Both reduce to closed forms in the partial moment phi(x) = E[(Y - x)^+] that
# gpd.R and gev.R already provide, so neither needs quadrature.
# ---------------------------------------------------------------------------

# --- (A) symmetric ---------------------------------------------------------
#
# For Y > t the weight is p and psi = max(Y - t, c); for Y < t the weight is
# (1-p) and -psi = max(t - Y, c). Splitting each at the knot,
#
#   E[max(Y-t, c) 1{Y>t}] = c (F(t+c) - F(t)) + phi(t+c) + c S(t+c)
#                         = c S(t) + phi(t+c)
#   E[max(t-Y, c) 1{Y<t}] = c (F(t) - F(t-c)) + phiL(t-c) + c F(t-c)
#                         = c F(t) + phiL(t-c)
#
# using F(t+c) + S(t+c) = 1, and phiL(x) = E[(x-Y)^+] = x - m + phi(x). So the
# identification equation is
#
#   p [c S(t) + phi(t+c)] = (1-p) [c F(t) + (t - c) - m + phi(t-c)].
#
# The model mean m survives, so (A) does NOT remove the mean anchoring.
gpd_invhuber_g <- function(t, p, mu, sigma, xi, cc) {
  m  <- gpd_mean(mu, sigma, xi)
  St <- sgpd(t, mu, sigma, xi)
  U  <- cc * St + gpd_partial_moment(t + cc, mu, sigma, xi)
  L  <- cc * (1 - St) + (t - cc) - m + gpd_partial_moment(t - cc, mu, sigma, xi)
  p * U - (1 - p) * L
}

# --- (B) one-sided ---------------------------------------------------------
#
#   p E[(Y-t)^+] = (1-p) E[min((t-Y)^+, c)]
#
# and E[min((t-Y)^+, c)] = phiL(t) - phiL(t-c) = c + phi(t) - phi(t-c), because
# the two mean terms cancel between phiL(t) and phiL(t-c). So
#
#   p phi(t) = (1-p) [c + phi(t) - phi(t-c)]
#
# is free of m entirely: capping the downside removes the mean anchoring
# structurally, rather than by substituting a noisy sample mean for it.
#
# Limits. As c -> Inf, phi(t-c) -> m - (t-c), the bracket -> phiL(t), and the
# expectile is recovered. As c -> 0 the bracket -> c F(t) -> 0, forcing
# phi(t) -> 0 and t to the upper endpoint: the family is degenerate at that end
# and does not reach the quantile. c is therefore a one-sided dial away from the
# expectile, not an interpolation between two usable estimators.
gpd_onesided_g <- function(t, p, mu, sigma, xi, cc) {
  ph <- gpd_partial_moment(t, mu, sigma, xi)
  p * ph - (1 - p) * (cc + ph - gpd_partial_moment(t - cc, mu, sigma, xi))
}

# --- shared solver ---------------------------------------------------------
# Both g are strictly decreasing in t (each is p*U(t) - (1-p)*L(t) with U
# decreasing and L increasing), so bisection from a bracket anchored on the
# quantile is safe and needs no derivative.
solve_decreasing <- function(g, lo, hi, tol = 1e-11, maxit = 200) {
  glo <- g(lo); ghi <- g(hi)
  for (it in 1:60) {                       # widen until it straddles
    bad <- glo <= 0
    if (!any(bad, na.rm = TRUE)) break
    span <- pmax(hi - lo, 1e-8)
    lo[bad] <- lo[bad] - span[bad]; glo <- g(lo)
  }
  for (it in 1:60) {
    bad <- ghi >= 0
    if (!any(bad, na.rm = TRUE)) break
    span <- pmax(hi - lo, 1e-8)
    hi[bad] <- hi[bad] + span[bad]; ghi <- g(hi)
  }
  for (it in seq_len(maxit)) {
    mid <- (lo + hi) / 2
    gm <- g(mid)
    up <- is.finite(gm) & gm > 0
    lo[up] <- mid[up]; hi[!up] <- mid[!up]
    if (max(hi - lo, na.rm = TRUE) < tol * max(1, max(abs(mid), na.rm = TRUE))) break
  }
  (lo + hi) / 2
}

gpd_invhuber <- function(p, mu, sigma, xi, cc) {
  if (!is.finite(gpd_mean(mu, sigma, xi))) return(rep(NA_real_, length(p)))
  q <- qgpd(p, mu, sigma, xi)
  s0 <- sigma + xi * (q - mu)
  solve_decreasing(function(t) gpd_invhuber_g(t, p, mu, sigma, xi, cc),
                   lo = q - s0, hi = q + s0)
}

# Safeguarded Newton. g is strictly decreasing and smooth away from the knot,
# and its derivative is closed form -- phi'(x) = -S(x), so
#
#   g'(t) = -p S(t) + (1-p) [S(t) - S(t-c)]
#
# which turns a ~40-step bisection into a ~6-step Newton. Bisection brackets are
# kept and used whenever a Newton step would leave them, so it cannot diverge.
gpd_onesided <- function(p, mu, sigma, xi, cc, tol = 1e-10, maxit = 40) {
  if (!is.finite(gpd_mean(mu, sigma, xi))) return(rep(NA_real_, length(p)))
  q <- qgpd(p, mu, sigma, xi)
  s0 <- sigma + xi * (q - mu)
  g <- function(t) gpd_onesided_g(t, p, mu, sigma, xi, cc)
  lo <- q - s0; hi <- q + s0
  glo <- g(lo); ghi <- g(hi)
  for (it in 1:60) {                          # widen until it straddles
    bad <- !is.finite(glo) | glo <= 0
    if (!any(bad)) break
    lo[bad] <- lo[bad] - pmax(hi - lo, 1e-8)[bad]; glo <- g(lo)
  }
  for (it in 1:60) {
    bad <- !is.finite(ghi) | ghi >= 0
    if (!any(bad)) break
    hi[bad] <- hi[bad] + pmax(hi - lo, 1e-8)[bad]; ghi <- g(hi)
  }
  # Newton converges in about five steps, but a handful of near-degenerate
  # levels need many more. Iterating the whole vector until the slowest one
  # settles costs every level the slowest one's budget, so converged elements
  # are frozen and dropped from the active set instead.
  x <- (lo + hi) / 2
  act <- rep(TRUE, length(x))
  for (it in seq_len(maxit)) {
    if (!any(act)) break
    i <- which(act)
    xi_ <- x[i]; pi_ <- p[i]
    gx <- gpd_onesided_g(xi_, pi_, mu, sigma, xi, cc)
    up <- is.finite(gx) & gx > 0
    lo[i[up]] <- xi_[up]; hi[i[!up]] <- xi_[!up]
    Sx <- sgpd(xi_, mu, sigma, xi)
    dg <- -pi_ * Sx + (1 - pi_) * (Sx - sgpd(xi_ - cc, mu, sigma, xi))
    step <- ifelse(is.finite(dg) & dg < 0, xi_ - gx / dg, (lo[i] + hi[i]) / 2)
    outside <- !is.finite(step) | step <= lo[i] | step >= hi[i]
    step[outside] <- ((lo[i] + hi[i]) / 2)[outside]
    done <- abs(step - xi_) <= tol * (1 + abs(xi_))
    x[i] <- step
    act[i[done]] <- FALSE
  }
  x
}

# --- generic versions, for an arbitrary distribution ------------------------
# Both equations are written in the partial moment phi(x) = E[(Y-x)^+] alone
# (family A additionally needs F and the mean), so they apply to any truth whose
# phi is available -- in particular the contaminated DGP, via
# partial_moment_true(). This is what lets the mean-anchoring claim be checked
# on the truth rather than only on the fitted family.
invhuber_general <- function(p, cc, phi, Fcdf, m, lo, hi, tol = 1e-11) {
  g <- function(t) {
    St <- 1 - Fcdf(t)
    p * (cc * St + phi(t + cc)) -
      (1 - p) * (cc * (1 - St) + (t - cc) - m + phi(t - cc))
  }
  solve_decreasing(g, lo = rep(lo, length(p)), hi = rep(hi, length(p)), tol = tol)
}

onesided_general <- function(p, cc, phi, lo, hi, tol = 1e-11) {
  g <- function(t) {
    ph <- phi(t)
    p * ph - (1 - p) * (cc + ph - phi(t - cc))
  }
  solve_decreasing(g, lo = rep(lo, length(p)), hi = rep(hi, length(p)), tol = tol)
}

# --- composite fitting ------------------------------------------------------
# The loss whose derivative is psi(u) = max(u, -c) is the "reverse Huber":
#
#   rho_c(u) = u^2                 for u >= -c      (quadratic on the upside)
#            = -2 c u - c^2        for u <  -c      (linear far below)
#
# continuous and C1 at u = -c. Weighting it by |p - I(y < t)| and integrating
# against w(p) gives the composite criterion. Splitting the sample at t and at
# t - c makes every piece a partial sum of y and y^2, so the whole grid is one
# pass over two cumulative sums -- the same trick the elastile fitter uses.
#
# c is fixed from the data before optimising (c = k * IQR) rather than tied to
# the fitted scale, so rho does not move with theta and the ordinary
# M-estimation theory applies unchanged. k is the tuning parameter.
gpd_onesided_loss <- function(tv, p, wts, ys, C1, C2, S1, S2, n, cc) {
  j  <- findInterval(tv, ys)                  # #{y <= t}
  jc <- findInterval(tv - cc, ys)             # #{y <= t - c}
  c1 <- C1[j + 1]; c2 <- C2[j + 1]
  d1 <- C1[jc + 1]; d2 <- C2[jc + 1]
  # y > t, weight p, quadratic
  hi  <- (S2 - c2) - 2 * tv * (S1 - c1) + (n - j) * tv^2
  # t - c <= y <= t, weight (1-p), quadratic
  mid <- (c2 - d2) - 2 * tv * (c1 - d1) + (j - jc) * tv^2
  # y < t - c, weight (1-p), linear
  low <- 2 * cc * (jc * tv - d1) - jc * cc^2
  sum(wts * (p * hi + (1 - p) * (mid + low)))
}

fit_gpd_onesided <- function(y, grid, w_fun, k) {
  cc <- k * max(1e-8, IQR(y))
  p <- grid$p; wts <- grid$w_quad * w_fun(p)
  keep <- wts > 0; p <- p[keep]; wts <- wts[keep]
  ys <- sort(y); n <- length(ys)
  C1 <- c(0, cumsum(ys)); C2 <- c(0, cumsum(ys^2)); S1 <- C1[n+1]; S2 <- C2[n+1]
  obj <- function(par) {
    mu <- par[1]; sigma <- exp(par[2]); xi <- par[3]
    if (!is.finite(mu) || !is.finite(sigma) || xi < XI_LO || xi > XI_HI) return(BIG)
    tv <- gpd_onesided(p, mu, sigma, xi, cc)
    if (any(!is.finite(tv))) return(BIG)
    gpd_onesided_loss(tv, p, wts, ys, C1, C2, S1, S2, n, cc)
  }
  start <- c(min(y) - 0.05 * diff(range(y)), log(max(1e-6, sd(y))), 0.15)
  o <- try(optim(start, obj, method = "Nelder-Mead",
                 control = list(maxit = 3000, reltol = 1e-12)), silent = TRUE)
  if (inherits(o, "try-error") || o$value >= BIG) return(rep(NA_real_, 3))
  o <- optim(o$par, obj, method = "Nelder-Mead",
             control = list(maxit = 3000, reltol = 1e-12))
  if (o$value >= BIG) return(rep(NA_real_, 3))
  c(o$par[1], exp(o$par[2]), o$par[3])
}

# ---------------------------------------------------------------------------
# Lp-quantiles of a GPD (Chen 1996; a special case of Breckling & Chambers'
# M-quantiles, and the object Daouia, Girard & Stupfler study in the extremes
# literature).
#
# The Lp-quantile at level p with exponent a minimises
#   E[ |p - I(Y < t)| |Y - t|^a ],
# so it solves the identification equation
#   p U_a(t) = (1 - p) L_a(t),
#   U_a(t) = E[(X - t)^{a-1} 1{X > t}],   L_a(t) = E[(t - X)^{a-1} 1{X < t}].
# a = 1 gives the quantile, a = 2 the expectile.
#
# For the GPD the UPPER piece is closed form. Substituting u = 1 + xi(x-mu)/sigma
# and then u = u_t / v turns it into a complete Beta integral:
#
#   U_a(t) = (sigma/xi)^{a-1} (1/xi) u_t^{a-1-1/xi} B(a, 1/xi - a + 1),
#
# finite exactly when 1/xi > a - 1. The LOWER piece runs over the bounded
# interval [mu, t] and has no such form, but the substitution x = t - (t-mu) v^2
# removes the endpoint singularity in (t-x)^{a-1} and leaves a smooth integrand,
# so a fixed Gauss-Legendre rule is accurate.
#
# Note the moment condition this family carries: the LOSS E|Y - t|^a is finite
# iff xi < 1/a, which recovers xi < 1 for the quantile (a = 1) and xi < 1/2 for
# the expectile (a = 2) as the two ends of one family.
# ---------------------------------------------------------------------------

gpd_upper_frac <- function(t, mu, sigma, xi, a) {
  # E[(X - t)^{a-1} 1{X > t}]
  stopifnot(xi > 0, 1 / xi > a - 1)
  out <- numeric(length(t))
  below <- t < mu
  ut <- ifelse(below, 1, 1 + xi * (t - mu) / sigma)
  ok <- ut > 0
  cf <- (sigma / xi)^(a - 1) * (1 / xi) * beta(a, 1 / xi - a + 1)
  out[ok] <- cf * ut[ok]^(a - 1 - 1 / xi)
  # Below the support the whole distribution lies above t; the integral picks up
  # the shifted contribution, which has no closed form, so fall back to quadrature.
  if (any(below)) out[below] <- gpd_frac_numeric(t[below], mu, sigma, xi, a, side = "upper")
  out
}

gpd_lower_frac <- function(t, mu, sigma, xi, a, n_gl = 48) {
  # E[(t - X)^{a-1} 1{X < t}], integrated in PROBABILITY space:
  #
  #   L(t) = int_0^{F(t)} (t - Q(u))^{a-1} du,   u = F(t) (1 - v^2),
  #
  # so L(t) = int_0^1 (t - Q(F(t)(1-v^2)))^{a-1} * 2 F(t) v dv. Integrating over u
  # rather than x keeps the mass evenly spread however far t sits from mu: the
  # earlier x-space substitution put a fixed number of nodes on [mu, t], so once
  # t >> mu (heavy xi) almost none of them landed where the density lives, and it
  # lost several digits. The v^2 substitution flattens the u -> F(t) endpoint,
  # where the integrand vanishes like v^{2a-1}. Vectorised over t as one
  # (n_gl x length(t)) matrix.
  gl <- gauss_legendre_01(n_gl)
  out <- numeric(length(t))
  ok <- t > mu
  if (!any(ok)) return(out)
  tt <- t[ok]
  Ft <- 1 - sgpd(tt, mu, sigma, xi)                   # length m
  v <- gl$x                                           # length n_gl
  U <- outer(1 - v^2, Ft)                             # u = F(t)(1 - v^2)
  Q <- qgpd(as.vector(U), mu, sigma, xi)
  dim(Q) <- dim(U)
  D <- pmax(rep(tt, each = n_gl) - Q, 0)
  W <- D^(a - 1) * 2 * rep(Ft, each = n_gl) * v
  out[ok] <- as.vector(gl$w %*% W)
  out
}

gpd_frac_numeric <- function(t, mu, sigma, xi, a, side = "upper", n_gl = 40) {
  gl <- gauss_legendre_01(n_gl)
  vapply(t, function(tt) {
    # integrate over probability, x = Q(u), u in (0, 1)
    s <- gl$x; u <- 1 - s^3
    x <- qgpd(u, mu, sigma, xi)
    val <- if (side == "upper") pmax(x - tt, 0)^(a - 1) else pmax(tt - x, 0)^(a - 1)
    sum(gl$w * val * 3 * s^2)
  }, numeric(1))
}

# The Lp-quantile function, by regula falsi (Illinois) on the identification
# equation, bracketed by the quantile and the expectile.
gpd_lp_quantile <- function(p, mu, sigma, xi, a, tol = 1e-10, maxit = 40) {
  if (abs(a - 1) < 1e-12) return(qgpd(p, mu, sigma, xi))
  if (abs(a - 2) < 1e-12) return(gpd_expectile(p, mu, sigma, xi))
  # Bracket. The expectile is only a hint, and it is undefined for xi >= 1 --
  # precisely the region where the L^a-quantile still exists (it needs only
  # xi < 1/(a-1)), so fall back to a scale-based bracket around the quantile.
  q <- qgpd(p, mu, sigma, xi)
  e <- if (xi < 1) gpd_expectile(p, mu, sigma, xi) else rep(NA_real_, length(q))
  e[!is.finite(e)] <- q[!is.finite(e)]
  s0 <- sigma + xi * (q - mu)                 # GPD scale at q; positive on the support
  lo <- pmin(q, e) - 1e-3 * s0; hi <- pmax(q, e) + 1e-3 * s0
  g <- function(t) p * gpd_upper_frac(t, mu, sigma, xi, a) -
                   (1 - p) * gpd_lower_frac(t, mu, sigma, xi, a)
  glo <- g(lo); ghi <- g(hi)
  # Widen if the bracket does not straddle (can happen at extreme levels).
  for (it in 1:40) {
    bad <- glo * ghi > 0
    if (!any(bad)) break
    span <- pmax(hi - lo, 1e-3 * s0)
    lo[bad] <- pmax(mu + 1e-12, lo[bad] - span[bad]); hi[bad] <- hi[bad] + span[bad]
    glo <- g(lo); ghi <- g(hi)
  }
  x <- lo; side <- integer(length(p))
  for (it in seq_len(maxit)) {
    x <- (lo * ghi - hi * glo) / (ghi - glo)
    bad <- !is.finite(x) | x <= lo | x >= hi
    x[bad] <- ((lo + hi) / 2)[bad]
    gx <- g(x)
    if (max(abs(gx)) < tol * max(1, max(abs(glo), abs(ghi)))) break
    up <- gx * ghi > 0                          # root in the lower half
    hi[up] <- x[up]; ghi[up] <- gx[up]
    glo[up & side == 1] <- glo[up & side == 1] / 2
    side[up] <- 1L
    lo[!up] <- x[!up]; glo[!up] <- gx[!up]
    ghi[!up & side == -1] <- ghi[!up & side == -1] / 2
    side[!up] <- -1L
    if (max(hi - lo) < tol * max(1, max(abs(x)))) break
  }
  x
}

# --- composite Lp-quantile estimator for the GPD ------------------------------
#
# Minimises  sum_i int w(p) |p - I(y_i < T_a(p|theta))| |y_i - T_a(p|theta)|^a dp
# with T_a the model's Lp-quantile function. Unlike the L1 and L2 cases the
# per-level sum cannot be reduced to cumulative sums, because |y - t|^a with
# non-integer a does not expand in powers of t; it is formed directly as an
# n-by-K matrix, which at these sizes costs little next to the Lp-quantile solve.
fit_gpd_lp_composite <- function(y, grid, w_fun, a) {
  p <- grid$p; wts <- grid$w_quad * w_fun(p)
  keep <- wts > 0; p <- p[keep]; wts <- wts[keep]
  n <- length(y)
  P <- matrix(p, nrow = n, ncol = length(p), byrow = TRUE)
  obj <- function(par) {
    mu <- par[1]; sigma <- exp(par[2]); xi <- par[3]
    if (!is.finite(mu) || !is.finite(sigma) || xi <= 1e-4 || xi > XI_HI) return(BIG)
    if (1 / xi <= a - 1) return(BIG)         # upper fractional moment must exist
    tv <- try(gpd_lp_quantile(p, mu, sigma, xi, a), silent = TRUE)
    if (inherits(tv, "try-error") || any(!is.finite(tv))) return(BIG)
    D <- outer(y, tv, "-")
    sum(wts * colSums(ifelse(D >= 0, P, 1 - P) * abs(D)^a))
  }
  start <- c(min(y) - 0.05 * diff(range(y)), log(max(1e-6, sd(y))), 0.15)
  o <- try(optim(start, obj, method = "Nelder-Mead",
                 control = list(maxit = 2000, reltol = 1e-10)), silent = TRUE)
  if (inherits(o, "try-error") || o$value >= BIG) return(rep(NA_real_, 3))
  o <- optim(o$par, obj, method = "Nelder-Mead", control = list(maxit = 2000, reltol = 1e-10))
  if (o$value >= BIG) return(rep(NA_real_, 3))
  c(o$par[1], exp(o$par[2]), o$par[3])
}

# ---------------------------------------------------------------------------
# The smooth graft evaluated directly, for speed.
#
# distplyr's eval_survival recomputes the correction integral from the anchor on
# every call, which is fine interactively and prohibitive inside a Monte Carlo
# (about 1.3 s per replicate for a return-period grid). Here the integral is
# accumulated once per replicate.
#
# The construction is the same one:
#   M(x) = (1-w) Sbody(x) + w Stail(x),
#   J(x) = int_{x0}^{x} w'(s) [Stail(s) - Sbody(s)] / M(s) ds,
#   S(x) = M(x) exp(-J(x)).
#
# Two structural facts make it cheap. The body is the empirical distribution, so
# Sbody is a step function and the integrand is smooth BETWEEN order statistics
# -- panel-by-panel Gauss-Legendre is then exact to quadrature precision. And
# w' vanishes outside the data range, so J is zero below the sample minimum and
# frozen above the maximum, where the graft is the fitted tail times a constant.
#
# Checked against distplyr's own eval_quantile in scripts/98-validate.R.
# ---------------------------------------------------------------------------

# Builds the graft's survival function and the frozen constant. Exposed so the
# construction itself can be checked against distplyr, not only the inversion.
graft_fast_build <- function(y, theta, w_fun, family = c("gev", "gpd"), n_gl = 8,
                             w_deriv = NULL) {
  family <- match.arg(family)
  if (any(is.na(theta))) return(NULL)
  ys <- sort(y); n <- length(ys)
  if (family == "gev") {
    Stail <- function(x) gev_survival(x, theta[1], theta[2], theta[3])
    Qtail <- function(p) qgev(p, theta[1], theta[2], theta[3])
    lo_tail <- gev_endpoints(theta[1], theta[2], theta[3])[1]
  } else {
    Stail <- function(x) sgpd(x, theta[1], theta[2], theta[3])
    Qtail <- function(p) qgpd(p, theta[1], theta[2], theta[3])
    lo_tail <- theta[1]
  }
  # Weight transported through the interpolated empirical cdf,
  #   w_x(x) = w(Fhat(x)),   w_x'(x) = w'(Fhat(x)) * Fhat'(x).
  # Fhat is piecewise linear with slope 1 / (n * gap) on each inter-order-statistic
  # panel, so with an analytic w' the transported derivative is exact -- and since
  # the quadrature panels are exactly those intervals, no panel ever straddles the
  # kinks. A central difference on x, by contrast, straddles them and is what makes
  # a numerically differentiated weight disagree in the fourth decimal.
  pp <- (seq_len(n) - 0.5) / n
  Fhat <- stats::approxfun(ys, pp, rule = 2)
  wx <- function(x) w_fun(Fhat(x))
  gaps <- diff(ys)
  if (is.null(w_deriv)) {
    hh <- 1e-6
    w_deriv <- function(p) (w_fun(pmin(1, p + hh)) - w_fun(pmax(0, p - hh))) /
                           (pmin(1, p + hh) - pmax(0, p - hh))
  }
  fhat <- function(x) {
    out <- numeric(length(x))
    inside <- x > ys[1] & x < ys[n]
    if (any(inside)) {
      i <- findInterval(x[inside], ys, all.inside = TRUE)
      out[inside] <- 1 / (n * gaps[i])
    }
    out
  }
  wpx <- function(x) w_deriv(Fhat(x)) * fhat(x)
  Sbody <- function(x) (n - findInterval(x, ys)) / n     # right-continuous
  M <- function(x) { w <- wx(x); (1 - w) * Sbody(x) + w * Stail(x) }
  integrand <- function(x) {
    m <- M(x)
    out <- wpx(x) * (Stail(x) - Sbody(x)) / m
    out[!is.finite(out)] <- 0
    out
  }
  # Accumulate J over the panels between order statistics.
  gl <- gauss_legendre_01(n_gl)
  Jknot <- numeric(n)                                    # J at each ys[i]
  if (n > 1) {
    lo <- ys[-n]; hi <- ys[-1]; wd <- hi - lo
    nodes <- outer(gl$x, wd) + rep(lo, each = n_gl)      # n_gl x (n-1)
    vals <- matrix(integrand(as.vector(nodes)), nrow = n_gl)
    panel <- as.vector(gl$w %*% vals) * wd
    Jknot[-1] <- cumsum(panel)
  }
  Jfun <- function(x) {                                  # J at arbitrary x
    x <- pmin(pmax(x, ys[1]), ys[n])
    i <- pmax(1, pmin(n - 1, findInterval(x, ys, all.inside = TRUE)))
    a <- ys[i]; b <- x
    wd <- b - a
    if (all(wd <= 0)) return(Jknot[i])
    nodes <- outer(gl$x, wd) + rep(a, each = n_gl)
    vals <- matrix(integrand(as.vector(nodes)), nrow = n_gl)
    Jknot[i] + as.vector(gl$w %*% vals) * wd
  }
  Sfun <- function(x) M(x) * exp(-Jfun(x))

  list(S = Sfun, cst = exp(-Jknot[n]), s_hi = Stail(ys[n]),
       Qtail = Qtail, ys = ys, Sknot = Sfun(ys))
}

graft_fast_return_levels <- function(y, theta, w_fun, ex, family = c("gev", "gpd"),
                                     n_gl = 8, w_deriv = NULL) {
  b <- graft_fast_build(y, theta, w_fun, family, n_gl, w_deriv)
  if (is.null(b)) return(rep(NA_real_, length(ex)))
  ys <- b$ys; n <- length(ys); Sfun <- b$S; Qtail <- b$Qtail
  cst <- b$cst; s_hi <- b$s_hi; Sknot <- b$Sknot
  out <- numeric(length(ex))
  for (k in seq_along(ex)) {
    e <- ex[k]
    if (is.finite(s_hi) && s_hi > 0 && e <= cst * s_hi) {
      out[k] <- Qtail(1 - e / cst)                       # tail regime, analytic
      next
    }
    i <- which(Sknot <= e)[1]
    if (is.na(i)) { out[k] <- ys[n]; next }
    if (i == 1) { out[k] <- ys[1]; next }
    # Either the atom at ys[i] straddles e, or the crossing is inside the panel.
    left <- Sfun(ys[i] - 1e-12)
    if (left <= e) {
      f <- function(x) Sfun(x) - e
      out[k] <- tryCatch(stats::uniroot(f, c(ys[i - 1], ys[i] - 1e-12), tol = 1e-10)$root,
                         error = function(...) ys[i])
    } else {
      out[k] <- ys[i]                                    # the jump crosses e
    }
  }
  attr(out, "c") <- cst
  out
}

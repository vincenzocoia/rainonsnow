# ---------------------------------------------------------------------------
# The data-generating process: a Gaussian-contaminated GEV.
#
# The truth is the distribution of  max(X, Z)  with X ~ GEV and Z ~ Normal
# independent, so its cdf is the product
#
#   F(x) = F_GEV(x; mu, sigma, xi) * Phi((x - a) / b).
#
# The normal sits in the middle of the distribution, so it dominates the body
# and vanishes in the tail: as x grows, Phi -> 1 and F -> F_GEV exactly. A GEV
# is therefore the correct tail model but a misspecified body model, which is
# the situation the composite estimators are meant to exploit.
# ---------------------------------------------------------------------------

# Default truth. Chosen in scripts/00-design.R; see that script for the
# diagnostics behind the normal's location and scale.
TRUTH <- list(mu = 0, sigma = 1, xi = 0.2, a = 1.5, b = 0.8)

p_true <- function(x, tr = TRUTH) {
  pgev(x, tr$mu, tr$sigma, tr$xi) * pnorm(x, tr$a, tr$b)
}

r_true <- function(n, tr = TRUTH) {
  # Exact: the product of cdfs is the cdf of the maximum of independent draws.
  pmax(qgev(runif(n), tr$mu, tr$sigma, tr$xi), rnorm(n, tr$a, tr$b))
}

q_true <- function(p, tr = TRUTH) {
  lo <- gev_endpoints(tr$mu, tr$sigma, tr$xi)[1]
  vapply(p, function(pp) {
    if (pp <= 0) return(lo)
    if (pp >= 1) return(Inf)
    # The true quantile is at least the GEV quantile (F_true <= F_GEV).
    hi <- qgev(pp, tr$mu, tr$sigma, tr$xi)
    hi <- max(hi, tr$a) + 1
    while (p_true(hi, tr) < pp) hi <- hi * 2 + 1
    uniroot(function(x) p_true(x, tr) - pp, c(lo + 1e-12, hi),
            tol = 1e-12)$root
  }, numeric(1))
}

# How much of the Gaussian contamination is left at x: the probability that the
# normal component still exceeds x. This is what decides where the weight
# function should switch on.
contamination <- function(x, tr = TRUTH) pnorm(x, tr$a, tr$b, lower.tail = FALSE)

# --- expectiles of the true (contaminated) distribution ---------------------
#
# Worth computing explicitly, because expectiles behave differently from
# quantiles here. The p-quantile of the truth equals the p-quantile of its GEV
# component as soon as the normal's cdf reaches 1: quantiles are tail-local.
# Expectiles are not. The p-expectile solves
#
#   k phi(x) + m - x = 0,   k = (2p-1)/(1-p),  phi(x) = E[(Y - x)^+],
#
# and while phi(x) for large x sees only the tail, the mean m sees the whole
# distribution -- including the contaminated body. So the truth's expectiles
# differ from its GEV component's expectiles at EVERY level, by an amount that
# only dies out as p -> 1. That is a structural limit on how unbiased a
# tail-weighted expectile estimator can be.

mean_true <- function(tr = TRUTH) {
  gl <- gauss_legendre_01(400)
  # substitution u = 1 - s^3 concentrates nodes near u = 1, where Q is steep
  s <- gl$x; u <- 1 - s^3
  sum(gl$w * q_true(u, tr) * 3 * s^2)
}

partial_moment_true <- function(x, tr = TRUTH) {
  # Split at X0, beyond which the normal component is numerically exhausted and
  # the truth's survival function equals the GEV's.
  X0 <- tr$a + 12 * tr$b
  vapply(x, function(xx) {
    tail_part <- gev_partial_moment(max(xx, X0), tr$mu, tr$sigma, tr$xi)
    if (xx >= X0) return(tail_part)
    body_part <- integrate(function(t) 1 - p_true(t, tr), xx, X0,
                           rel.tol = 1e-11, subdivisions = 2000)$value
    body_part + tail_part
  }, numeric(1))
}

expectile_true <- function(p, tr = TRUTH, m = NULL) {
  if (is.null(m)) m <- mean_true(tr)
  lo_end <- gev_endpoints(tr$mu, tr$sigma, tr$xi)[1]
  vapply(p, function(pp) {
    if (pp == 0.5) return(m)
    k <- (2 * pp - 1) / (1 - pp)
    g <- function(x) k * partial_moment_true(x, tr) + m - x
    lo <- if (pp > 0.5) m else lo_end
    hi <- if (pp > 0.5) m + 1 else m
    while (pp > 0.5 && g(hi) > 0) hi <- m + (hi - m) * 4
    uniroot(g, c(lo, hi), tol = 1e-10)$root
  }, numeric(1))
}

# ---------------------------------------------------------------------------
# Fitting a GPD, two ways.
#
# The competitor an astute practitioner would actually use: peaks over
# threshold. Pick a threshold, fit a two-parameter GPD to the exceedances by
# maximum likelihood or L-moments, and HARD graft it onto the empirical body at
# that threshold. Return levels come from the standard POT formula.
#
# The composite alternative: treat the GPD as a three-parameter family for the
# whole distribution and fit it by the composite criterion with a weight that
# rises from p = 0, then SMOOTH graft it onto the empirical body using that same
# weight as the handover. Here, unlike the GEV case, tying the two weights
# together is the natural choice: the composite weight already starts at zero on
# the body, which is exactly where the empirical distribution should be trusted.
# ---------------------------------------------------------------------------

# --- peaks over threshold ----------------------------------------------------

gpd_exceed_nll <- function(par, z) {
  sigma <- exp(par[1]); xi <- par[2]
  if (xi < XI_LO || xi > XI_HI) return(BIG)
  u <- 1 + xi * z / sigma
  if (any(u <= 0)) return(BIG)
  length(z) * log(sigma) + (1 + 1 / xi) * sum(log(u))
}

fit_pot_mle <- function(y, v) {
  u <- as.numeric(stats::quantile(y, v, type = 1))
  z <- y[y > u] - u
  if (length(z) < 5) return(c(u = u, sigma = NA, xi = NA, zeta = NA))
  o <- try(optim(c(log(mean(z)), 0.1), gpd_exceed_nll, z = z, method = "Nelder-Mead",
                 control = list(maxit = 2000, reltol = 1e-12)), silent = TRUE)
  if (inherits(o, "try-error") || o$value >= BIG) return(c(u = u, sigma = NA, xi = NA, zeta = NA))
  c(u = u, sigma = exp(o$par[1]), xi = o$par[2], zeta = length(z) / length(y))
}

fit_pot_lmom <- function(y, v) {
  u <- as.numeric(stats::quantile(y, v, type = 1))
  z <- sort(y[y > u] - u); m <- length(z)
  if (m < 5) return(c(u = u, sigma = NA, xi = NA, zeta = NA))
  # First two sample L-moments of the exceedances.
  b0 <- mean(z); b1 <- sum((seq_len(m) - 1) / (m - 1) * z) / m
  l1 <- b0; l2 <- 2 * b1 - b0
  # For a GPD with location 0: l1 = sigma/(1-xi), l2 = sigma/((1-xi)(2-xi)),
  # so xi = 2 - l1/l2 and sigma = l1 (1 - xi).
  xi <- 2 - l1 / l2
  sigma <- l1 * (1 - xi)
  if (!is.finite(xi) || !is.finite(sigma) || sigma <= 0 || xi > XI_HI || xi < XI_LO)
    return(c(u = u, sigma = NA, xi = NA, zeta = NA))
  c(u = u, sigma = sigma, xi = xi, zeta = m / length(y))
}

# Return levels of a hard graft: empirical below the threshold, fitted GPD
# above, joined so that P(X > u) is the empirical exceedance rate zeta.
pot_return_levels <- function(y, fit, ex) {
  u <- fit[["u"]]; sigma <- fit[["sigma"]]; xi <- fit[["xi"]]; zeta <- fit[["zeta"]]
  out <- numeric(length(ex))
  emp <- distionary::dst_empirical(y)
  for (i in seq_along(ex)) {
    e <- ex[i]
    if (is.na(sigma) || e >= zeta) {
      out[i] <- as.numeric(distionary::eval_quantile(emp, at = 1 - e))
    } else if (abs(xi) < 1e-12) {
      out[i] <- u + sigma * log(zeta / e)
    } else {
      out[i] <- u + sigma * ((e / zeta)^(-xi) - 1) / xi
    }
  }
  out
}

# --- composite estimation of a three-parameter GPD ---------------------------

make_gpd_composite_objective <- function(y, grid, w_fun, type, anchored = FALSE) {
  p <- grid$p; wts <- grid$w_quad * w_fun(p)
  keep <- wts > 0; p <- p[keep]; wts <- wts[keep]
  ys <- sort(y); n <- length(ys)
  C1 <- c(0, cumsum(ys)); C2 <- c(0, cumsum(ys^2))
  S1 <- C1[n + 1]; S2 <- C2[n + 1]
  m_hat <- mean(y)
  wp <- wts * p; wq <- wts * (1 - p)
  function(par) {
    mu <- par[1]; sigma <- exp(par[2]); xi <- par[3]
    if (!is.finite(mu) || !is.finite(sigma) || xi < XI_LO || xi > XI_HI) return(BIG)
    tv <- if (type == "quantile") qgpd(p, mu, sigma, xi)
          else gpd_expectile(p, mu, sigma, xi, if (anchored) m_hat else NULL)
    if (any(!is.finite(tv))) return(BIG)
    j <- findInterval(tv, ys); c1 <- C1[j + 1]; c2 <- C2[j + 1]
    if (type == "quantile") {
      sum(wts * (p * (S1 - n * tv) - (c1 - j * tv)))
    } else {
      lo <- c2 - 2 * tv * c1 + j * tv^2
      hi <- (S2 - c2) - 2 * tv * (S1 - c1) + (n - j) * tv^2
      sum(wq * lo + wp * hi)
    }
  }
}

fit_gpd_composite <- function(y, grid, w_fun, type, anchored = FALSE) {
  # Start from a moment-ish guess: location just below the sample minimum.
  start <- c(min(y) - 0.05 * diff(range(y)), log(max(1e-6, sd(y))), 0.15)
  obj <- make_gpd_composite_objective(y, grid, w_fun, type, anchored)
  o <- try(optim(start, obj, method = "Nelder-Mead",
                 control = list(maxit = 3000, reltol = 1e-12)), silent = TRUE)
  if (inherits(o, "try-error") || o$value >= BIG) return(rep(NA_real_, 3))
  o <- optim(o$par, obj, method = "Nelder-Mead", control = list(maxit = 3000, reltol = 1e-12))
  if (o$value >= BIG) return(rep(NA_real_, 3))
  c(o$par[1], exp(o$par[2]), o$par[3])
}

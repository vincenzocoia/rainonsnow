# ---------------------------------------------------------------------------
# Four ways of fitting a GEV.
#
#   fit_mle       maximum likelihood
#   fit_lmom      L-moments (Hosking)
#   fit_cqe       composite QUANTILE estimator   (L1; the thesis estimator)
#   fit_cee       composite EXPECTILE estimator  (L2; the idea under test)
#
# All return c(mu, sigma, xi), or NA's on failure.
# ---------------------------------------------------------------------------

XI_LO <- -0.45      # keep the fitted family in a regime where the mean exists
XI_HI <-  0.90      # and the numerics are stable
BIG   <- 1e12

# --- maximum likelihood ----------------------------------------------------

gev_nll <- function(par, y) {
  mu <- par[1]; sigma <- exp(par[2]); xi <- par[3]
  if (xi < XI_LO || xi > XI_HI) return(BIG)
  if (abs(xi) < 1e-8) {
    u <- (y - mu) / sigma
    return(length(y) * log(sigma) + sum(u + exp(-u)))
  }
  t <- 1 + xi * (y - mu) / sigma
  if (any(t <= 0)) return(BIG)          # data outside the fitted support
  z <- t^(-1 / xi)
  length(y) * log(sigma) + (1 + 1 / xi) * sum(log(t)) + sum(z)
}

fit_mle <- function(y, start = NULL) {
  if (is.null(start)) start <- fit_lmom(y)
  if (any(is.na(start))) start <- c(mean(y), sd(y), 0.1)
  par0 <- c(start[1], log(max(start[2], 1e-6)), min(max(start[3], -0.3), 0.6))
  # Make sure the starting value is feasible; MLE needs all data in support.
  if (gev_nll(par0, y) >= BIG) par0 <- c(mean(y), log(sd(y)), 0.05)
  o <- try(optim(par0, gev_nll, y = y, method = "Nelder-Mead",
                 control = list(maxit = 2000, reltol = 1e-12)), silent = TRUE)
  if (inherits(o, "try-error") || o$value >= BIG) return(rep(NA_real_, 3))
  o <- optim(o$par, gev_nll, y = y, method = "Nelder-Mead",
             control = list(maxit = 2000, reltol = 1e-12))
  c(o$par[1], exp(o$par[2]), o$par[3])
}

# --- L-moments (Hosking 1990) ----------------------------------------------

sample_lmoments <- function(y) {
  x <- sort(y); n <- length(x); i <- seq_len(n)
  b0 <- mean(x)
  b1 <- sum((i - 1) / (n - 1) * x) / n
  b2 <- sum((i - 1) * (i - 2) / ((n - 1) * (n - 2)) * x) / n
  c(l1 = b0, l2 = 2 * b1 - b0, l3 = 6 * b2 - 6 * b1 + b0)
}

# Hosking's GEV is F(x) = exp(-(1 - k (x - loc)/scale)^(1/k)), so his k is the
# negative of the shape used here: xi = -k.
lmom_to_gev <- function(l1, l2, l3) {
  t3 <- l3 / l2
  cc <- 2 / (3 + t3) - log(2) / log(3)
  k <- 7.8590 * cc + 2.9554 * cc^2
  if (!is.finite(k) || abs(k) < 1e-9) {          # Gumbel limit
    scale <- l2 / log(2)
    loc <- l1 - scale * 0.577215664901533
    return(c(loc, scale, 0))
  }
  gk <- gamma(1 + k)
  scale <- l2 * k / ((1 - 2^(-k)) * gk)
  loc <- l1 - scale * (1 - gk) / k
  c(loc, scale, -k)
}

fit_lmom <- function(y) {
  lm <- sample_lmoments(y)
  out <- try(lmom_to_gev(lm[1], lm[2], lm[3]), silent = TRUE)
  if (inherits(out, "try-error") || any(!is.finite(out)) || out[2] <= 0)
    return(rep(NA_real_, 3))
  unname(out)
}

# --- composite estimators ---------------------------------------------------
#
# Both minimise  sum_i \int_0^1 w(p) rho_p(y_i, T(p | theta)) dp  where T is the
# quantile function (L1) or the expectile function (L2), and
#
#   quantile:  rho_p(y, q) = (p - I(y < q)) (y - q)              [pinball]
#   expectile: rho_p(y, e) = |p - I(y < e)| (y - e)^2            [asymmetric LS]
#
# Each rho is the canonical consistent scoring function for its functional, so
# the integrated loss is minimised in expectation by a theta whose quantile
# (resp. expectile) function matches the truth wherever w > 0.
#
# The integral is a fixed Gauss-Legendre rule (see R/quadrature.R), so each
# estimator is a well-defined M-estimator, identical across replicates.

# The two losses, written so that the sum over observations at each level is a
# pair of table lookups rather than an n x K matrix.
#
# With y sorted and cumulative sums C1[j] = sum_{i<=j} y_i, C2[j] = sum_{i<=j} y_i^2,
# and j(e) = #{i : y_i < e}, the per-level sums are
#
#   sum_{y_i <  e} (y_i - e)^2 = C2[j] - 2 e C1[j] + j e^2
#   sum_{y_i >= e} (y_i - e)^2 = (S2 - C2[j]) - 2 e (S1 - C1[j]) + (n - j) e^2
#   sum_i (p - I(y_i < e))(y_i - e) = p (S1 - n e) - (C1[j] - j e)
#
# so the loss is exact, not approximated -- only the arithmetic changes.
make_composite_objective <- function(y, grid, w_fun, type) {
  p <- grid$p
  wts <- grid$w_quad * w_fun(p)
  keep <- wts > 0
  p <- p[keep]; wts <- wts[keep]
  ys <- sort(y); n <- length(ys)
  C1 <- c(0, cumsum(ys)); C2 <- c(0, cumsum(ys^2))
  S1 <- C1[n + 1]; S2 <- C2[n + 1]
  wp <- wts * p; wq <- wts * (1 - p)
  function(par) {
    mu <- par[1]; sigma <- exp(par[2]); xi <- par[3]
    if (!is.finite(mu) || !is.finite(sigma) || xi < XI_LO || xi > XI_HI)
      return(BIG)
    tv <- if (type == "quantile") qgev(p, mu, sigma, xi)
          else gev_expectile(p, mu, sigma, xi)
    if (any(!is.finite(tv))) return(BIG)
    j <- findInterval(tv, ys)          # #{i : y_i <= tv}; ties are measure zero
    c1 <- C1[j + 1]; c2 <- C2[j + 1]
    if (type == "quantile") {
      sum(wts * (p * (S1 - n * tv) - (c1 - j * tv)))
    } else {
      lo <- c2 - 2 * tv * c1 + j * tv^2                       # sum below
      hi <- (S2 - c2) - 2 * tv * (S1 - c1) + (n - j) * tv^2   # sum at/above
      sum(wq * lo + wp * hi)
    }
  }
}

fit_composite <- function(y, grid, w_fun, type, start = NULL) {
  if (is.null(start)) start <- fit_lmom(y)
  if (any(is.na(start))) start <- c(mean(y), sd(y), 0.1)
  obj <- make_composite_objective(y, grid, w_fun, type)
  par0 <- c(start[1], log(max(start[2], 1e-6)), min(max(start[3], XI_LO + .05), 0.6))
  o <- try(optim(par0, obj, method = "Nelder-Mead",
                 control = list(maxit = 3000, reltol = 1e-12)), silent = TRUE)
  if (inherits(o, "try-error") || o$value >= BIG) return(rep(NA_real_, 3))
  # One restart: Nelder-Mead can stall on a ridge, and the tail-focused loss
  # surface is poorly conditioned in (mu, sigma).
  o <- optim(o$par, obj, method = "Nelder-Mead",
             control = list(maxit = 3000, reltol = 1e-12))
  if (o$value >= BIG) return(rep(NA_real_, 3))
  c(o$par[1], exp(o$par[2]), o$par[3])
}

fit_cqe <- function(y, grid, w_fun, start = NULL)
  fit_composite(y, grid, w_fun, "quantile", start)
fit_cee <- function(y, grid, w_fun, start = NULL)
  fit_composite(y, grid, w_fun, "expectile", start)

# --- the alpha-elastile composite estimator ---------------------------------
#
# Minimises  sum_i int w(p) rho^alpha_p(y_i, T_alpha(p | theta)) dp  with
#
#   rho^alpha_p(y, t) = (alpha/s) |p - I(y<t)| (y-t)^2
#                     + (1 - alpha) (p - I(y<t)) (y - t)
#
# and T_alpha the model's alpha-elastile function. The sum of two strictly
# consistent scoring functions evaluated at the same point is strictly
# consistent for the minimiser of the sum, so this is a proper criterion
# targeting the alpha-elastile, with alpha = 0 and 1 the two estimators already
# studied.
#
# The scale s makes the two terms commensurable -- without it alpha = 1/2 is not
# a midpoint, since the L2 term has units of y^2 and the L1 term units of y. It
# is fixed per dataset at the ratio of the two integrated losses evaluated at
# the preliminary L-moment fit, so that at alpha = 1/2 the two contribute
# equally there. Being a constant during the optimisation, it does not disturb
# the properness of the criterion; it only chooses the units in which alpha
# interpolates.

elastile_scale <- function(y, grid, w_fun, start) {
  p <- grid$p; wts <- grid$w_quad * w_fun(p)
  keep <- wts > 0; p <- p[keep]; wts <- wts[keep]
  ys <- sort(y); n <- length(ys)
  C1 <- c(0, cumsum(ys)); C2 <- c(0, cumsum(ys^2))
  S1 <- C1[n + 1]; S2 <- C2[n + 1]
  piece <- function(tv) {
    j <- findInterval(tv, ys); c1 <- C1[j + 1]; c2 <- C2[j + 1]
    lo <- c2 - 2 * tv * c1 + j * tv^2
    hi <- (S2 - c2) - 2 * tv * (S1 - c1) + (n - j) * tv^2
    c(L2 = sum(wts * ((1 - p) * lo + p * hi)),
      L1 = sum(wts * (p * (S1 - n * tv) - (c1 - j * tv))))
  }
  e <- gev_expectile(p, start[1], start[2], start[3])
  q <- qgev(p, start[1], start[2], start[3])
  s <- piece(e)[["L2"]] / piece(q)[["L1"]]
  if (!is.finite(s) || s <= 0) s <- max(1e-8, sd(y))
  s
}

make_elastile_objective <- function(y, grid, w_fun, alpha, s) {
  p <- grid$p; wts <- grid$w_quad * w_fun(p)
  keep <- wts > 0; p <- p[keep]; wts <- wts[keep]
  ys <- sort(y); n <- length(ys)
  C1 <- c(0, cumsum(ys)); C2 <- c(0, cumsum(ys^2))
  S1 <- C1[n + 1]; S2 <- C2[n + 1]
  function(par) {
    mu <- par[1]; sigma <- exp(par[2]); xi <- par[3]
    if (!is.finite(mu) || !is.finite(sigma) || xi < XI_LO || xi > XI_HI) return(BIG)
    tv <- gev_elastile(p, mu, sigma, xi, alpha, s)
    if (any(!is.finite(tv))) return(BIG)
    j <- findInterval(tv, ys); c1 <- C1[j + 1]; c2 <- C2[j + 1]
    lo <- c2 - 2 * tv * c1 + j * tv^2
    hi <- (S2 - c2) - 2 * tv * (S1 - c1) + (n - j) * tv^2
    L2 <- sum(wts * ((1 - p) * lo + p * hi))
    L1 <- sum(wts * (p * (S1 - n * tv) - (c1 - j * tv)))
    (alpha / s) * L2 + (1 - alpha) * L1
  }
}

fit_elastile <- function(y, grid, w_fun, alpha, start = NULL) {
  if (is.null(start)) start <- fit_lmom(y)
  if (any(is.na(start))) start <- c(mean(y), sd(y), 0.1)
  s <- elastile_scale(y, grid, w_fun, start)
  obj <- make_elastile_objective(y, grid, w_fun, alpha, s)
  par0 <- c(start[1], log(max(start[2], 1e-6)), min(max(start[3], XI_LO + .05), 0.6))
  o <- try(optim(par0, obj, method = "Nelder-Mead",
                 control = list(maxit = 3000, reltol = 1e-12)), silent = TRUE)
  if (inherits(o, "try-error") || o$value >= BIG) return(rep(NA_real_, 3))
  o <- optim(o$par, obj, method = "Nelder-Mead",
             control = list(maxit = 3000, reltol = 1e-12))
  if (o$value >= BIG) return(rep(NA_real_, 3))
  c(o$par[1], exp(o$par[2]), o$par[3])
}

# --- mean-anchored composite expectile estimator ----------------------------
#
# Identical to fit_cee() except that the model's expectile function is computed
# with the SAMPLE mean in place of the model's own mean. The expectile
# identification equation splits into a tail-local part (the partial moment phi,
# where a GEV fitted to a contaminated sample is right) and a global part (the
# mean, where it is not). Taking the second from the data removes the residual
# that no weight function can.
#
# The price: this is no longer a strictly consistent scoring function for the
# model's expectile. It is a two-step estimating equation targeting a modified
# functional -- the mean is profiled out nonparametrically -- and needs its own
# justification. Fisher consistency survives: if the model is correct over the
# weighted region and the sample mean is consistent for the true mean, the true
# parameter still solves it.

make_anchored_objective <- function(y, grid, w_fun) {
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
    tv <- gev_expectile_anchored(p, mu, sigma, xi, m_hat)
    if (any(!is.finite(tv))) return(BIG)
    j <- findInterval(tv, ys); c1 <- C1[j + 1]; c2 <- C2[j + 1]
    lo <- c2 - 2 * tv * c1 + j * tv^2
    hi <- (S2 - c2) - 2 * tv * (S1 - c1) + (n - j) * tv^2
    sum(wq * lo + wp * hi)
  }
}

fit_cee_anchored <- function(y, grid, w_fun, start = NULL) {
  if (is.null(start)) start <- fit_lmom(y)
  if (any(is.na(start))) start <- c(mean(y), sd(y), 0.1)
  obj <- make_anchored_objective(y, grid, w_fun)
  par0 <- c(start[1], log(max(start[2], 1e-6)), min(max(start[3], XI_LO + .05), 0.6))
  o <- try(optim(par0, obj, method = "Nelder-Mead",
                 control = list(maxit = 3000, reltol = 1e-12)), silent = TRUE)
  if (inherits(o, "try-error") || o$value >= BIG) return(rep(NA_real_, 3))
  o <- optim(o$par, obj, method = "Nelder-Mead", control = list(maxit = 3000, reltol = 1e-12))
  if (o$value >= BIG) return(rep(NA_real_, 3))
  c(o$par[1], exp(o$par[2]), o$par[3])
}

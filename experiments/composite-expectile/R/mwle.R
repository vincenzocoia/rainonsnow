# ---------------------------------------------------------------------------
# Maximum weighted likelihood (Fung 2022, Insurance: Mathematics and Economics
# 107, 180-198), adapted from finite mixtures to a single extreme-value family.
#
# The estimator maximises  sum_i w(y_i) log f(y_i)  where f is the WEIGHT-TILTED
# density  f(y) = h(y) w(y) / C(theta),  C(theta) = int h(u;theta) w(u) du.
# Since log w(y_i) does not involve theta it drops, leaving
#
#     l(theta) = sum_i w(y_i) log h(y_i; theta)  -  (sum_i w(y_i)) log C(theta).
#
# The normaliser is what makes this different from a naive weighted likelihood:
# without it, tilting the weights towards the tail would simply rescale the
# objective. C(theta) is the expectation of the weight under the model.
#
# MAPPING THE WEIGHT. Fung's w is a function of the observation y; the composite
# criterion's weight is a function of the probability level p. They are put on
# the same footing by composing with the empirical distribution function,
# w(y) = w_p(Fhat(y)). Because Fhat is a step function this also makes C(theta)
# exact rather than quadrature: Fhat is j/n on [y_(j), y_(j+1)), so
#
#     C(theta) = sum_{j=0}^{n} w_p(j/n) [ H(y_(j+1); theta) - H(y_(j); theta) ]
#
# with y_(0) = -Inf and y_(n+1) = +Inf.
#
# Fung requires w bounded away from zero, so the composite weight is floored:
# w_p(p) -> eps + (1 - eps) w_p(p).
# ---------------------------------------------------------------------------

mwle_weights <- function(y, w_fun, eps = 1e-3) {
  n <- length(y)
  # Fhat at each observation, using the mid-rank convention for ties-free data
  Fh <- rank(y, ties.method = "max") / n
  eps + (1 - eps) * w_fun(Fh)
}

# log C(theta) computed exactly from the step structure of Fhat
mwle_logC <- function(ys, cdf_fun, w_fun, eps = 1e-3) {
  n <- length(ys)
  wj <- eps + (1 - eps) * w_fun(seq.int(0, n) / n)      # w on each flat piece
  Hy <- c(0, cdf_fun(ys), 1)                            # H at -Inf, data, +Inf
  sum(wj * diff(Hy))
}

# Generic fitter. `dens` and `cdf` take (x, par); `par` is on its natural scale.
fit_mwle <- function(y, w_fun, dens, cdf, start, lower_ok, eps = 1e-3) {
  ys <- sort(y)
  wi <- mwle_weights(y, w_fun, eps)
  sw <- sum(wi)
  nll <- function(par) {
    if (!lower_ok(par)) return(BIG)
    d <- dens(y, par)
    if (any(!is.finite(d)) || any(d <= 0)) return(BIG)
    C <- mwle_logC(ys, function(x) cdf(x, par), w_fun, eps)
    if (!is.finite(C) || C <= 0) return(BIG)
    v <- -(sum(wi * log(d)) - sw * log(C))
    if (!is.finite(v)) BIG else v
  }
  o <- try(optim(start, nll, method = "Nelder-Mead",
                 control = list(maxit = 3000, reltol = 1e-12)), silent = TRUE)
  if (inherits(o, "try-error") || o$value >= BIG) return(rep(NA_real_, length(start)))
  o <- optim(o$par, nll, method = "Nelder-Mead",
             control = list(maxit = 3000, reltol = 1e-12))
  if (o$value >= BIG) return(rep(NA_real_, length(start)))
  o$par
}

# --- GEV, three parameters (mu, log sigma, xi) ------------------------------
fit_mwle_gev <- function(y, w_fun, start = NULL, eps = 1e-3) {
  if (is.null(start)) start <- fit_lmom(y)
  if (any(is.na(start))) start <- c(mean(y), sd(y), 0.1)
  p0 <- c(start[1], log(max(start[2], 1e-6)), min(max(start[3], XI_LO + .05), 0.6))
  r <- fit_mwle(y, w_fun,
    dens = function(x, p) dgev(x, p[1], exp(p[2]), p[3]),
    cdf  = function(x, p) pgev(x, p[1], exp(p[2]), p[3]),
    start = p0,
    lower_ok = function(p) is.finite(p[1]) && is.finite(p[2]) &&
                            p[3] > XI_LO && p[3] < XI_HI, eps = eps)
  if (anyNA(r)) return(rep(NA_real_, 3))
  c(r[1], exp(r[2]), r[3])
}

# --- GPD with the threshold known at zero, two parameters -------------------
fit_mwle_gpd2 <- function(y, w_fun, start = NULL, eps = 1e-3) {
  if (is.null(start)) start <- gpd2_lmom(y)
  if (any(is.na(start))) start <- c(sd(y), 0.15)
  p0 <- c(log(max(start[1], 1e-6)), min(max(start[2], XI_LO + .05), 0.6))
  r <- fit_mwle(y, w_fun,
    dens = function(x, p) dgpd(x, 0, exp(p[1]), p[2]),
    cdf  = function(x, p) pgpd(x, 0, exp(p[1]), p[2]),
    start = p0,
    lower_ok = function(p) is.finite(p[1]) && p[2] > XI_LO && p[2] < XI_HI,
    eps = eps)
  if (anyNA(r)) return(rep(NA_real_, 2))
  c(exp(r[1]), r[2])
}

# --- GPD, three parameters (threshold estimated) ----------------------------
fit_mwle_gpd3 <- function(y, w_fun, eps = 1e-3) {
  st <- c(min(y) - 0.05 * diff(range(y)), log(max(1e-6, sd(y))), 0.15)
  r <- fit_mwle(y, w_fun,
    dens = function(x, p) dgpd(x, p[1], exp(p[2]), p[3]),
    cdf  = function(x, p) pgpd(x, p[1], exp(p[2]), p[3]),
    start = st,
    lower_ok = function(p) is.finite(p[1]) && p[1] < min(y) &&
                            is.finite(p[2]) && p[3] > XI_LO && p[3] < XI_HI,
    eps = eps)
  if (anyNA(r)) return(rep(NA_real_, 3))
  c(r[1], exp(r[2]), r[3])
}

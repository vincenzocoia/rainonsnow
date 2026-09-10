# ---------------------------------------------------------------------------
# Extremiles (Daouia, Gijbels & Stupfler).
#
# For tau >= 1/2 put r(tau) = log(1/2)/log(tau) and K_tau(t) = t^r, so the
# weight-generating function is J_tau(t) = K'_tau(t) = r t^(r-1). The extremile
# is the L-functional
#     xi_tau = int_0^1 J_tau(t) Q(t) dt = E[Y J_tau(F(Y))],
# equal to E[max(Y_1, ..., Y_r)] when r is an integer, and existing whenever
# E|Y| < infinity.
#
# Both extreme-value families give it in closed form.
#   GPD(0, s, xi):  Q(t) = s[(1-t)^-xi - 1]/xi, and r int_0^1 t^(r-1)(1-t)^-xi dt
#                   = Gamma(r+1)Gamma(1-xi)/Gamma(r+1-xi), so
#                   xi_tau = (s/xi)[Gamma(r+1)Gamma(1-xi)/Gamma(r+1-xi) - 1].
#   GEV(m, s, xi):  substituting u = -log t turns the integral into a gamma
#                   function directly, giving
#                   xi_tau = m + (s/xi)[r^xi Gamma(1-xi) - 1].
# Both need xi < 1. At tau = 1/2 (r = 1) each reduces to the family's mean.
# ---------------------------------------------------------------------------

extremile_r <- function(tau) log(0.5) / log(tau)      # tau >= 1/2

gpd_extremile <- function(tau, mu, sigma, xi) {
  if (xi >= 1) return(rep(NA_real_, length(tau)))
  r <- extremile_r(tau)
  if (abs(xi) < 1e-10) return(mu + sigma * (digamma(r + 1) - digamma(1)))
  mu + (sigma / xi) * (exp(lgamma(r + 1) + lgamma(1 - xi) - lgamma(r + 1 - xi)) - 1)
}

gev_extremile <- function(tau, mu, sigma, xi) {
  if (xi >= 1) return(rep(NA_real_, length(tau)))
  r <- extremile_r(tau)
  if (abs(xi) < 1e-10) return(mu + sigma * (log(r) + EULER))
  mu + (sigma / xi) * (r^xi * gamma(1 - xi) - 1)
}

# Empirical extremile: a linear combination of order statistics with weights
# K_tau(i/n) - K_tau((i-1)/n).
extremile_emp <- function(tau, y) {
  ys <- sort(y); n <- length(ys); i <- seq_len(n)
  vapply(tau, function(tt) {
    r <- extremile_r(tt)
    sum(((i / n)^r - ((i - 1) / n)^r) * ys)
  }, numeric(1))
}

# ---------------------------------------------------------------------------
# Composite extremile estimation.
#
# Extremiles are L-functionals, not M-quantiles: the weight J_tau(F(Y)) depends
# on the unknown F rather than on the residual, so there is no scoring function
# whose expected minimiser is the extremile, and the composite LOSS construction
# does not transfer. What does transfer is minimum-distance matching on the
# extremile function itself,
#
#     theta-hat = argmin  int w(tau) [ xi-hat_tau - xi_tau(theta) ]^2 dtau,
#
# with the empirical extremile a linear combination of order statistics and the
# model extremile in closed form. This is a tail-weighted probability-weighted-
# moment estimator, so it sits beside Hosking and Wallis rather than beside the
# M-estimators.
#
# THE UPPER LIMIT MATTERS, and it is the trap of section 7 in another guise. As
# tau -> 1, r(tau) -> infinity and xi_tau tends to the upper endpoint, while the
# empirical extremile saturates: it is an average of order statistics and cannot
# reach past the sample. Measured on GPD(0, 0.6, 0.25), the relative bias of the
# empirical extremile depends on r/n and does NOT vanish with n:
#
#     r/n     0.02    0.05    0.10    0.20    0.40    1.00
#     n=100  -0.001  -0.008  -0.025  -0.037  -0.071  -0.149
#     n=1000 -0.003  -0.009  -0.018  -0.028  -0.056  -0.129
#
# So a grid capped at r = n -- the first thing one would try -- puts the heaviest
# weight exactly where the empirical extremile is worst, and the estimator is not
# consistent. The cap must instead be an intermediate sequence: r_max = alpha n
# with alpha small and, for consistency, alpha -> 0 as n grows. This is the
# familiar k -> infinity, k/n -> 0 condition of extreme-value estimation, and it
# is the price of the extremile being an L-functional: it cannot reach as far
# into the tail as an empirical quantile can.
# ---------------------------------------------------------------------------

extremile_grid <- function(n, n_gl = 96, tau_lo = 0.5, alpha = 0.05) {
  tau_hi <- 0.5^(1 / (alpha * n))             # r(tau_hi) = alpha * n
  gl <- gauss_legendre_01(n_gl)
  list(tau = tau_lo + (tau_hi - tau_lo) * gl$x,
       w_quad = (tau_hi - tau_lo) * gl$w, tau_hi = tau_hi, r_hi = alpha * n)
}

fit_composite_extremile <- function(y, grid, w_fun, family = c("gev", "gpd2", "gpd"),
                                    start = NULL) {
  family <- match.arg(family)
  tau <- grid$tau; wt <- grid$w_quad * w_fun(tau)
  xh <- extremile_emp(tau, y)
  mod <- switch(family,
    gev  = function(p) gev_extremile(tau, p[1], exp(p[2]), p[3]),
    gpd  = function(p) gpd_extremile(tau, p[1], exp(p[2]), p[3]),
    gpd2 = function(p) gpd_extremile(tau, 0, exp(p[1]), p[2]))
  ok <- switch(family,
    gev  = function(p) is.finite(p[1]) && p[3] > XI_LO && p[3] < XI_HI,
    gpd  = function(p) is.finite(p[1]) && p[1] < min(y) && p[3] > XI_LO && p[3] < XI_HI,
    gpd2 = function(p) p[2] > XI_LO && p[2] < XI_HI)
  obj <- function(p) {
    if (!ok(p)) return(BIG)
    m <- mod(p)
    if (any(!is.finite(m))) return(BIG)
    v <- sum(wt * (xh - m)^2)
    if (is.finite(v)) v else BIG
  }
  if (is.null(start)) start <- switch(family,
    gev  = { s <- fit_lmom(y); if (anyNA(s)) c(mean(y), sd(y), 0.1) else s },
    gpd  = c(min(y) - 0.05 * diff(range(y)), max(1e-6, sd(y)), 0.15),
    gpd2 = { s <- gpd2_lmom(y); if (anyNA(s)) c(sd(y), 0.15) else s })
  p0 <- switch(family,
    gev  = c(start[1], log(max(start[2], 1e-6)), min(max(start[3], XI_LO + .05), 0.6)),
    gpd  = c(start[1], log(max(start[2], 1e-6)), min(max(start[3], XI_LO + .05), 0.6)),
    gpd2 = c(log(max(start[1], 1e-6)), min(max(start[2], XI_LO + .05), 0.6)))
  o <- try(optim(p0, obj, method = "Nelder-Mead",
                 control = list(maxit = 3000, reltol = 1e-12)), silent = TRUE)
  if (inherits(o, "try-error") || o$value >= BIG) return(rep(NA_real_, length(p0)))
  o <- optim(o$par, obj, method = "Nelder-Mead",
             control = list(maxit = 3000, reltol = 1e-12))
  if (o$value >= BIG) return(rep(NA_real_, length(p0)))
  switch(family,
    gev  = c(o$par[1], exp(o$par[2]), o$par[3]),
    gpd  = c(o$par[1], exp(o$par[2]), o$par[3]),
    gpd2 = c(exp(o$par[1]), o$par[2]))
}

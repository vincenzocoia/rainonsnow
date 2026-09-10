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

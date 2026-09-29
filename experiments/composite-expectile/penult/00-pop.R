# ---------------------------------------------------------------------------
# Populations, by their reciprocal hazard r(x) = S(x)/f(x).
#
# A GPD anchored at u has r(u+y) = sigma + xi y exactly: a straight line whose
# intercept is the scale and whose slope is the shape. So r is the natural
# coordinate for penultimate approximation -- r(x) is the penultimate scale at
# level x and r'(x) the penultimate shape, and the ultimate shape is lim r'.
#
# Both families below have r and r' in closed form, derived rather than
# differenced, because the whole analysis differentiates things built on them.
# ---------------------------------------------------------------------------

## ---- GEV(0, 1, xi) --------------------------------------------------------
# t(x) = (1+xi x)^(-1/xi);  S = 1-exp(-t);  f = t^(1+xi) exp(-t)
#   r(x)  = (exp(t) - 1) / t^(1+xi)
#   r'(x) = (1+xi)(exp(t) - 1)/t - exp(t)
# As t -> 0, r'(x) = xi - t(1-xi)/2 + O(t^2): r' rises to xi from BELOW.
pop_gev <- function(xi) {
  tfun <- if (abs(xi) < 1e-12) function(x) exp(-x) else
    function(x) (1 + xi * x)^(-1/xi)
  em1 <- function(t) expm1(t)                       # exp(t)-1, accurate at t~0
  list(
    name = sprintf("GEV xi=%.2f", xi), xi_inf = xi,
    S = function(x) -expm1(-tfun(x)),
    Q = function(p) if (abs(xi) < 1e-12) -log(-log(p)) else ((-log(p))^(-xi) - 1)/xi,
    ## quantile at a SURVIVAL level, so the far tail keeps its digits: at
    ## s = 1e-10, Q(1-s) loses them all to cancellation in log(p).
    Qs = function(s) { t <- -log1p(-s)
                       if (abs(xi) < 1e-12) -log(t) else (t^(-xi) - 1)/xi },
    r = function(x) { t <- tfun(x); em1(t) / t^(1 + xi) },
    rp = function(x) { t <- tfun(x); (1 + xi) * em1(t) / t - exp(t) }
  )
}

## ---- LP3: log X ~ Gamma(shape a, scale b) ---------------------------------
# m(t) = S_G(t)/f_G(t) is the Gamma's reciprocal hazard, and with t = log x
#   r(x)  = x m(t)
#   r'(x) = m(t) [1 + 1/b - (a-1)/t] - 1
# Expanding m(t) = b[1 + (a-1)b/t + O(t^-2)] gives r'(x) = b + b^2(a-1)/log x,
# which is the logarithmic (rho = 0) convergence this family is known for.
pop_lp3 <- function(a, b) {
  m <- function(t) exp(pgamma(t, shape = a, scale = b, lower.tail = FALSE, log.p = TRUE) -
                       dgamma(t, shape = a, scale = b, log = TRUE))
  list(
    name = sprintf("LP3 a=%.1f b=%.2f", a, b), xi_inf = b,
    S = function(x) pgamma(log(x), shape = a, scale = b, lower.tail = FALSE),
    Q = function(p) exp(qgamma(p, shape = a, scale = b)),
    Qs = function(s) exp(qgamma(s, shape = a, scale = b, lower.tail = FALSE)),
    r = function(x) x * m(log(x)),
    rp = function(x) { t <- log(x); m(t) * (1 + 1/b - (a - 1)/t) - 1 }
  )
}

## ---- a GPD, for unit tests: r must come back a straight line ---------------
pop_gpd <- function(sigma, xi, mu = 0) list(
  name = sprintf("GPD(%.2f,%.2f)", sigma, xi), xi_inf = xi,
  S = function(x) if (abs(xi) < 1e-12) exp(-(x-mu)/sigma) else (1 + xi*(x-mu)/sigma)^(-1/xi),
  Q = function(p) mu + sigma * ((1 - p)^(-xi) - 1)/xi,
  Qs = function(s) if (abs(xi) < 1e-12) mu - sigma*log(s) else mu + sigma*(s^(-xi) - 1)/xi,
  r = function(x) sigma + xi * (x - mu),
  rp = function(x) rep(xi, length(x))
)

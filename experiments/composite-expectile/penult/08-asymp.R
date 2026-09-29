# ---------------------------------------------------------------------------
# Asymptotic normality of the composite M-quantile estimators.
#
# THE STRUCTURE. The criterion is
#     S_n(theta) = int_0^1 w(p) (1/n) sum_i rho_p(Y_i - T_p(theta)) dp.
# Differentiating in theta -- the indicator contributes nothing, because rho is
# continuous at zero -- theta-hat solves (1/n) sum_i g(Y_i, theta-hat) = 0 with
#
#     g(y, theta) = int_0^1 w(p) h_p(y, theta) grad_theta T_p(theta) dp,
#     h_p(y, theta) = |p - I(y < T_p(theta))| psi(y - T_p(theta)).
#
# So it is a Z-estimator with a score that is itself an integral over levels,
# and the standard sandwich applies: sqrt(n)(theta-hat - theta*) -> N(0, A^-1 B A^-T)
# with theta* the root of E_F g = 0 -- exactly the pseudo-true parameter of
# section 17 -- and
#
#     A = -E[ d g / d theta' ]
#       = int w(p) [ c_p grad T_p grad T_p' - d_p Hess T_p ] dp
#     B = Var[ g(Y, theta*) ]
#       = int int w(p) w(q) kappa(p,q) grad T_p grad T_q' dp dq
#
#     c_p = E[ |p - I(Y < t_p)| psi'(Y - t_p) ]        (a local sensitivity)
#     d_p = E[ |p - I(Y < t_p)| psi(Y - t_p) ]         (the DISCREPANCY, = 0
#                                                       under correct
#                                                       specification)
#     kappa(p,q) = E[ h_p(Y) h_q(Y) ]                  (a covariance kernel
#                                                       across levels)
#
# Two readings worth having before any numbers. The Hessian term in A carries
# d_p, which is the same level-by-level discrepancy that defines theta* -- so it
# vanishes identically when the model is right, and A collapses to the outer
# product alone. And B is a DOUBLE integral: the estimator's variance is driven
# by how correlated the level-p and level-q criteria are, which is where the
# kernel width of section 21 shows up in the variance rather than the bias.
# ---------------------------------------------------------------------------
source("penult/04-setup.R")

## c_p, in closed form for each psi. psi' is a delta for the pinball loss, which
## is why the quantile row is a density rather than a probability.
c_of_p <- function(ex, kind, p, t, cc = NULL, alpha = NULL, s_el = NULL) {
  Ft <- ex$Fc(t); St <- 1 - Ft
  dens <- function(z) ex$Sc(z) / ex$rq(z)               # f_c = S_c / r
  switch(kind,
    quantile  = dens(t),
    expectile = p * St + (1 - p) * Ft,
    onesided  = p * St + (1 - p) * (Ft - ex$Fc(t - cc)),
    elastile  = (2 * alpha / s_el) * (p * St + (1 - p) * Ft) + (1 - alpha) * dens(t))
}

## kappa(p,q) = E[h_p h_q]. The pinball loss gets the closed form, because
## h_p(y) = p - I(y < t_p) is a step and quadrature of a step on a fixed grid is
## only O(de) accurate. Everything else is integrated: h_p is continuous there,
## and one matrix product gives the whole kernel at once.
kappa_matrix <- function(ex, kind, p, t, cc = NULL, alpha = NULL, s_el = NULL) {
  if (kind == "quantile") {
    Ft <- ex$Fc(t); n <- length(p)
    Fmin <- outer(t, t, function(a, b) ex$Fc(pmin(a, b)))
    return(outer(p, p) - outer(p, Ft) - outer(Ft, p) + Fmin)
  }
  y <- ex$y; a <- ex$e
  N <- length(a); wt <- rep(c(2, 4), length.out = N); wt[1] <- 1; wt[N] <- 1
  meas <- wt * ex$de / 3 * exp(-ex$e) * (ex$r0 / ex$r)     # the truth's measure
  psi <- switch(kind,
    expectile = function(v) v,
    onesided  = function(v) pmax(v, -cc),
    elastile  = function(v) (2 * alpha / s_el) * v + (1 - alpha) * sign(v))
  H <- vapply(seq_along(p), function(i) {
    lam <- ifelse(y > t[i], p[i], 1 - p[i])
    lam * psi(y - t[i]) }, numeric(N))                     # N x n_levels
  crossprod(H * meas, H)                                   # kappa, n x n
}

## the sandwich
sandwich_cov <- function(ex, kind, w_fun, grid, theta, cc = NULL, alpha = NULL,
                         s_el = NULL) {
  Tf <- mk_T(kind, cc, alpha, s_el); gf <- mk_g(ex, kind, cc, alpha, s_el)
  p <- grid$p; wq <- grid$w_quad * w_fun(p)
  k <- wq > 1e-14; p <- p[k]; wq <- wq[k]
  t <- Tf(p, theta)
  ## grad and Hessian of T in theta, by central differences. The step must sit
  ## well above the functional solvers' own tolerance (1e-14 here); the lesson
  ## from the target solver is that nesting differences too tightly turns the
  ## derivative into amplified solver noise.
  h <- 1e-4 * pmax(abs(theta), 1e-2)
  Tpm <- function(d) Tf(p, theta + d)
  G1 <- (Tpm(c(h[1], 0)) - Tpm(c(-h[1], 0))) / (2 * h[1])
  G2 <- (Tpm(c(0, h[2])) - Tpm(c(0, -h[2]))) / (2 * h[2])
  H11 <- (Tpm(c(h[1], 0)) - 2 * t + Tpm(c(-h[1], 0))) / h[1]^2
  H22 <- (Tpm(c(0, h[2])) - 2 * t + Tpm(c(0, -h[2]))) / h[2]^2
  H12 <- (Tpm(c(h[1], h[2])) - Tpm(c(h[1], -h[2])) -
          Tpm(c(-h[1], h[2])) + Tpm(c(-h[1], -h[2]))) / (4 * h[1] * h[2])
  cp <- c_of_p(ex, kind, p, t, cc, alpha, s_el)
  dp <- gf(t, p)
  A <- matrix(0, 2, 2)
  A[1,1] <- sum(wq * (cp * G1 * G1 - dp * H11))
  A[2,2] <- sum(wq * (cp * G2 * G2 - dp * H22))
  A[1,2] <- A[2,1] <- sum(wq * (cp * G1 * G2 - dp * H12))
  KM <- kappa_matrix(ex, kind, p, t, cc, alpha, s_el)
  V <- cbind(wq * G1, wq * G2)
  B <- crossprod(V, KM %*% V)
  Ai <- solve(A)
  list(A = A, B = B, V = Ai %*% B %*% t(Ai), dp_max = max(abs(dp)))
}

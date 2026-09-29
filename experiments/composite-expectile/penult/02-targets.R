# ---------------------------------------------------------------------------
# Population targets: the GPD each estimator converges to above a fixed u.
#
# Replace the empirical distribution in the criterion by the true exceedance
# law and optimise. Nothing here is Monte Carlo.
#
# For the composite family the target solves the estimating equation derived in
# the main report,
#     int w(p) g_p(T_p(theta)) dT_p/dtheta dp = 0,
#     g_p(t) = E_Fc[ |p - I(Z < t)| psi(Z - t) ],
# where g_p is the truth's own M-quantile identification function evaluated at
# the MODEL's level-p functional. It vanishes exactly when the two agree, so
# theta* sets a sensitivity-weighted average of level-by-level discrepancies to
# zero. That is what makes the estimator a line fit through r, and what the
# kernel below measures.
# ---------------------------------------------------------------------------

## ---- a damped Newton for the 2-d system, warm-startable --------------------
# Restarted Nelder-Mead (what the earlier pseudo-true script used) is far too
# slow to sit inside the kernel loop, which re-solves the target ~80 times per
# estimator. Warm-started Newton converges in 2-4 steps.
# NOTE ON STEP SIZES. G already contains a finite difference (dT/dtheta), so
# differencing G again is a nested difference: with an outer step near or below
# the inner one, the inner solver's own tolerance is amplified by 1/h_outer and
# the Jacobian comes out numerical noise -- in testing, the wrong sign. The
# outer step must sit two decades above the inner one, and the functionals must
# be solved tighter than either. h_in below is the inner step, set by the caller.
solve2 <- function(G, start, tol = 1e-11, maxit = 80, lo = c(1e-8, -0.9), hi = c(Inf, 0.95)) {
  th <- start
  for (it in seq_len(maxit)) {
    g0 <- G(th); if (any(!is.finite(g0))) return(list(par = rep(NA_real_, 2), resid = Inf))
    if (sqrt(sum(g0^2)) < tol) break
    h <- 1e-3 * pmax(abs(th), 1e-2)
    J <- cbind((G(th + c(h[1], 0)) - G(th - c(h[1], 0))) / (2 * h[1]),
               (G(th + c(0, h[2])) - G(th - c(0, h[2]))) / (2 * h[2]))
    st <- try(solve(J, -g0), silent = TRUE)
    if (inherits(st, "try-error")) return(list(par = rep(NA_real_, 2), resid = Inf))
    ok <- FALSE
    for (lam in c(1, 0.5, 0.25, 0.1, 0.03)) {          # backtrack on ||G||
      cand <- pmin(pmax(th + lam * st, lo), hi)
      gc_ <- G(cand)
      if (all(is.finite(gc_)) && sqrt(sum(gc_^2)) < sqrt(sum(g0^2))) { th <- cand; ok <- TRUE; break }
    }
    if (!ok) break
  }
  rz <- sqrt(sum(G(th)^2))
  list(par = th, resid = if (is.finite(rz)) rz else Inf)
}

## ---- expectations against the exceedance law ------------------------------
# Two integrators, and they are not interchangeable. The grid is uniform in the
# BASE level a; a perturbed law's own level e_new(a) is not, so an expectation
# under the truth carries the extra factor d e_new/d a = r0/r_new. Using the
# plain rule for it silently biases every perturbed target, which is exactly
# where a kernel would go wrong without announcing it.
mk_simp <- function(N, de) { wt <- rep(c(2, 4), length.out = N); wt[1] <- 1; wt[N] <- 1
                             wt * de / 3 }
mk_intT <- function(ex) {                    # E under the truth's exceedance law
  N <- length(ex$a); wt <- mk_simp(N, ex$de)
  jac <- exp(-ex$e) * (ex$r0 / ex$r)
  function(vals) sum(wt * jac * vals)
}
mk_intU <- function(ex) {                    # int_0^Inf exp(-b) (.) db, model side
  N <- length(ex$a); wt <- mk_simp(N, ex$de); eb <- exp(-ex$a)
  function(vals) sum(wt * eb * vals)
}

## ---- references ------------------------------------------------------------
target_mle <- function(ex, start = c(1, 0.2)) {             # KL projection
  I <- mk_intT(ex); y <- ex$y
  nll <- function(par) {
    s <- exp(par[1]); xi <- par[2]
    z <- 1 + xi * y / s
    if (any(z <= 0)) return(1e10)
    -I(-log(s) - (1/xi + 1) * log(z))
  }
  o <- optim(c(log(start[1]), start[2]), nll, control = list(reltol = 1e-14, maxit = 4000))
  o <- optim(o$par, nll, control = list(reltol = 1e-14, maxit = 4000))
  c(exp(o$par[1]), o$par[2])
}
target_lmom <- function(ex) {                               # L-moment matching
  I <- mk_intT(ex); y <- ex$y
  l1 <- I(y); l2 <- I(y * (1 - 2 * exp(-ex$e)))
  xi <- 2 - l1 / l2; c(l1 * (1 - xi), xi)
}

## ---- the composite family --------------------------------------------------
mk_g <- function(ex, kind, cc = NULL, alpha = NULL, s_el = NULL) {
  phiL <- mk_phiL(ex)
  switch(kind,
    quantile  = function(t, p) p - ex$Fc(t),
    expectile = function(t, p) p * ex$phic(t) - (1 - p) * phiL(t),
    elastile  = function(t, p) (2 * alpha / s_el) * (p * ex$phic(t) - (1 - p) * phiL(t)) +
                               (1 - alpha) * (p - ex$Fc(t)),
    onesided  = function(t, p) p * ex$phic(t) - (1 - p) * (phiL(t) - phiL(t - cc)))
}
# The functionals are solved to 1e-14 rather than their defaults: they sit two
# finite differences below the Jacobian, so their tolerance sets its noise floor.
mk_T <- function(kind, cc = NULL, alpha = NULL, s_el = NULL) switch(kind,
  quantile  = function(p, th) qgpd(p, 0, th[1], th[2]),
  expectile = function(p, th) gpd_expectile(p, 0, th[1], th[2], tol = 1e-14, maxit = 200),
  elastile  = function(p, th) gpd_elastile(p, 0, th[1], th[2], alpha, s_el),
  onesided  = function(p, th) gpd_onesided(p, 0, th[1], th[2], cc, tol = 1e-14, maxit = 200))

target_composite <- function(ex, kind, w_fun, grid, start = c(1, 0.2),
                             cc = NULL, alpha = NULL, s_el = NULL) {
  gf <- mk_g(ex, kind, cc, alpha, s_el); Tf <- mk_T(kind, cc, alpha, s_el)
  p <- grid$p; wq <- grid$w_quad * w_fun(p)
  k <- wq > 1e-14; p <- p[k]; wq <- wq[k]; wq <- wq / sum(wq)
  G <- function(th) {
    if (th[1] <= 0 || th[2] <= -0.85 || th[2] >= 0.95) return(c(NA, NA))
    tv <- Tf(p, th); if (any(!is.finite(tv))) return(c(NA, NA))
    gv <- gf(tv, p)
    h <- 1e-5 * pmax(abs(th), 1e-2)
    d1 <- (Tf(p, th + c(h[1], 0)) - Tf(p, th - c(h[1], 0))) / (2 * h[1])
    d2 <- (Tf(p, th + c(0, h[2])) - Tf(p, th - c(0, h[2]))) / (2 * h[2])
    c(sum(wq * gv * d1) / th[1], sum(wq * gv * d2))       # scaled for conditioning
  }
  ## The one-sided equation has a spurious root near the shape boundary -- it
  ## cost a day of the finite-sample study -- so a result that lands there, or
  ## fails to drive |G| down, is retried from a spread of starts.
  bad <- function(r) isTRUE(!all(is.finite(r$par))) || !isTRUE(r$resid <= 1e-8) ||
                     isTRUE(r$par[2] < -0.6) || isTRUE(r$par[2] > 0.9)
  r <- solve2(G, start)
  if (bad(r)) for (st in list(c(start[1], 0.05), c(start[1], 0.35), c(start[1], 0.6),
                              c(start[1] * 2, 0.2), c(start[1] / 2, 0.2))) {
    r2 <- solve2(G, st)
    if (!bad(r2)) { r <- r2; break }
    if (is.finite(r2$resid) && r2$resid < r$resid) r <- r2
  }
  r$par
}

## ---- Fung's weighted likelihood, population version ------------------------
# omega acts on the value, w on the level; bridged through the truth's own Fc,
# omega(y) = w(Fc(y)), which is how the finite-sample study mapped them.
target_mwle <- function(ex, w_fun, start = c(1, 0.2)) {
  I <- mk_intT(ex); IU <- mk_intU(ex); y <- ex$y; a <- ex$e
  wv <- w_fun(1 - exp(-a)); Wbar <- I(wv)          # = int_0^1 w(v) dv
  nll <- function(par) {
    s <- exp(par[1]); xi <- par[2]
    z <- 1 + xi * y / s; if (any(z <= 0)) return(1e10)
    lg <- -log(s) - (1/xi + 1) * log(z)
    ## C(theta) = int g(y) omega(y) dy = int_0^1 omega(Q_theta(q)) dq
    b <- ex$a                                  # model-side level, not the truth's
    yq <- if (abs(xi) < 1e-9) s * b else s * expm1(xi * b) / xi
    Cth <- IU(w_fun(ex$Fc(yq)))
    -(I(wv * lg) - Wbar * log(max(Cth, 1e-300)))
  }
  o <- optim(c(log(start[1]), start[2]), nll, control = list(reltol = 1e-14, maxit = 4000))
  o <- optim(o$par, nll, control = list(reltol = 1e-14, maxit = 4000))
  c(exp(o$par[1]), o$par[2])
}

## ---- composite extremile ---------------------------------------------------
# An L-functional, so there is no loss to integrate; the target is the
# minimum-distance match on the extremile function, over a grid whose reach
# r_max stands in for the intermediate sequence alpha*n.
target_extremile <- function(ex, w_fun, r_max = 50, n_gl = 64, start = c(1, 0.2)) {
  I <- mk_intT(ex); y <- ex$y; a <- ex$e
  gl <- gauss_legendre_01(n_gl)
  rr <- exp(log(1) + gl$x * (log(r_max) - log(1)))            # log-spaced reach
  tau <- 0.5^(1 / rr)
  wt <- gl$w * (log(r_max)) * rr * w_fun(tau); wt <- wt / sum(wt)
  emp <- vapply(rr, function(R) I(R * (-expm1(-a))^(R - 1) * y), numeric(1))
  obj <- function(par) {
    s <- exp(par[1]); xi <- par[2]; if (xi >= 0.95) return(1e10)
    mod <- gpd_extremile(tau, 0, s, xi)
    if (any(!is.finite(mod))) return(1e10)
    sum(wt * (emp - mod)^2)
  }
  o <- optim(c(log(start[1]), start[2]), obj, control = list(reltol = 1e-14, maxit = 4000))
  o <- optim(o$par, obj, control = list(reltol = 1e-14, maxit = 4000))
  c(exp(o$par[1]), o$par[2])
}

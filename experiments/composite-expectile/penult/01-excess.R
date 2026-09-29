# ---------------------------------------------------------------------------
# The exceedance law above u, built from r alone.
#
# Everything an estimator's population criterion needs -- F, the quantile
# function, both partial moments -- follows from r(u+y), because
#
#     e(y) = -log Sc(y) = int_0^y dt / r(u+t).
#
# The first version integrated that as an ODE on a uniform e-grid. This one
# fixes a BASE grid once, at uniform e for the unperturbed population, and then
# changes variable to the base level a: with t = y(a) and dt = r0(u+y(a)) da,
#
#     e_new(y(A)) = int_0^A [r0/r_new](u + y(a)) da
#     phi_new(y(A)) = int_A^Inf exp(-e_new(a)) r0(u + y(a)) da
#
# -- uniform-grid integrals of smooth integrands that sit near 1, vectorised
# instead of stepped. Two gains beyond speed: the unperturbed run reproduces
# e = a exactly, and every perturbed run in a kernel sweep uses the identical
# grid, so the quadrature error cancels in the difference rather than adding to
# it.
# ---------------------------------------------------------------------------

## cumulative integral of f on a uniform grid, O(h^4) (Adams-Moulton form).
## Plain cumulative trapezoid is O(h^2) and not enough: the targets are wanted
## to 1e-6 and are differenced again for the kernel.
cumsimp <- function(f, h) {
  n <- length(f); out <- numeric(n)
  inc <- h/12 * (5*f[1:(n-2)] + 8*f[2:(n-1)] - f[3:n])
  lastinc <- h/12 * (-f[n-2] + 8*f[n-1] + 5*f[n])
  out[-1] <- cumsum(c(inc, lastinc))
  out
}
revcumsimp <- function(f, h) rev(cumsimp(rev(f), h))     # int from the right

## The base grid: y placed at uniform exceedance level for the true population.
make_base <- function(P, u, Emax = 20, N = 6001) {
  a <- seq(0, Emax, length.out = N); de <- a[2] - a[1]
  Su <- P$S(u)
  y <- P$Qs(Su * exp(-a)) - u; y[1] <- 0
  list(P = P, u = u, Su = Su, a = a, de = de, y = y, r0 = P$r(u + y),
       N = N, Emax = Emax)
}

## rfun = NULL gives the unperturbed law
excess_on <- function(base, rfun = NULL) {
  y <- base$y; a <- base$a; de <- base$de; N <- base$N; u <- base$u
  rn <- if (is.null(rfun)) base$r0 else rfun(u + y)
  e <- cumsimp(base$r0 / rn, de)
  Sn <- exp(-e)
  g <- Sn * base$r0
  ## beyond the grid, continue as the GPD (r, r') define at the top edge
  xi_hi <- (rn[N] - rn[N-1]) / (y[N] - y[N-1])
  phi <- revcumsimp(g, de) + Sn[N] * rn[N] / max(1 - xi_hi, 1e-3)
  EZ <- phi[1]
  l <- log1p(y)
  sp_e   <- splinefun(l, e, method = "monoH.FC")
  sp_phi <- splinefun(l, log(pmax(phi, 1e-300)), method = "monoH.FC")
  sp_y   <- splinefun(e, l, method = "monoH.FC")
  y_hi <- y[N]; sig_hi <- rn[N]; phi_hi <- phi[N]; S_hi <- Sn[N]
  ## the naive (1 + xi d/sigma)^(-1/xi) silently returns 1 as xi -> 0
  gpd_pow <- function(d, pw) if (abs(xi_hi) < 1e-8) exp(-pw * d / sig_hi) else
    pmax(1 + xi_hi * d / sig_hi, 0)^(-pw / xi_hi)
  Scfun <- function(z) {
    o <- numeric(length(z)); hi <- z > y_hi
    o[!hi] <- exp(-sp_e(log1p(pmax(z[!hi], 0))))
    if (any(hi)) o[hi] <- S_hi * gpd_pow(z[hi] - y_hi, 1)
    pmin(pmax(o, 0), 1) }
  list(u = u, EZ = EZ, e = e, y = y, r = rn, r0 = base$r0, de = de, a = a, ymax = y_hi,
       Sc = Scfun, Fc = function(z) 1 - Scfun(z),
       phic = function(z) {
         o <- numeric(length(z)); hi <- z > y_hi; lo <- z <= 0; k <- !hi & !lo
         o[k] <- exp(sp_phi(log1p(z[k]))); o[lo] <- EZ
         if (any(hi)) o[hi] <- phi_hi * gpd_pow(z[hi] - y_hi, 1 - xi_hi)
         o },
       Qc = function(p) expm1(sp_y(-log1p(-p))),
       ## r and the density on the exceedance law, for the sandwich's c_p term
       rq = function(z) approx(y, rn, pmin(pmax(z, 0), y_hi), rule = 2)$y,
       fc = function(z) Scfun(z) / approx(y, rn, pmin(pmax(z, 0), y_hi), rule = 2)$y,
       e_of_y = function(z) sp_e(log1p(pmax(z, 0))))
}
# E[(z-Z)^+] = z - E[Z] + E[(Z-z)^+] above the support's lower end, and is
# identically 0 at or below it. Returning z for z <= 0 sends the one-sided
# estimating equation to a spurious boundary root.
mk_phiL <- function(ex) function(z) ifelse(z <= 0, 0, z - ex$EZ + ex$phic(pmax(z, 0)))

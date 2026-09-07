# ---------------------------------------------------------------------------
# 29. Pseudo-true limits derived, not simulated.
#
# scripts/28 estimated the pseudo-true parameters from 2e5 exceedances. That was
# not precise enough to separate the three composite estimators, whose fitted
# shapes differed by about the same amount as the seed-to-seed noise. This gets
# them exactly, from the population estimating equations.
#
# THE DERIVATION. The composite criterion is
#     S(theta) = int w(p) E_F[ |p - I(Y<T(p|theta))| rho(Y - T(p|theta)) ] dp.
# Differentiating in theta -- the indicator contributes nothing, because rho is
# continuous at zero -- the pseudo-true theta* solves
#
#     int w(p) g_p(T(p|theta)) dT(p|theta)/dtheta dp = 0,          (*)
#
#     g_p(t) = E_F[ |p - I(Y < t)| psi(Y - t) ].
#
# g_p(t) is exactly the M-quantile identification function of the TRUTH,
# evaluated at the MODEL's level-p functional. It vanishes iff the model's
# functional equals the truth's at that level. So (*) says: the pseudo-true
# theta sets a weighted average of the level-by-level discrepancies to zero,
# weighted by w(p) and by the sensitivity dT/dtheta. Two immediate consequences:
# the target depends on the weight w, and it depends on psi only through which
# functional is being matched.
#
# For the truth's conditional law above u, writing Fc, phic, phiLc for its cdf
# and upper/lower partial moments, the identification functions are
#   quantile   g = p - Fc(t)
#   expectile  g = p phic(t) - (1-p) phiLc(t)
#   elastile   g = (2a/s)[p phic - (1-p) phiLc] + (1-a)[p - Fc(t)]
#   one-sided  g = p phic(t) - (1-p)[phiLc(t) - phiLc(t-c)]
#
# References are obtained exactly too: the MLE limit minimises Kullback-Leibler
# divergence, the L-moment limit matches the first two L-moments.
#
# Output: out/pseudotrue-exact.txt
# ---------------------------------------------------------------------------
source("R/setup.R"); source("R/config.R")
source("R/gpd.R"); source("R/gpd_estimators.R"); source("R/invhuber.R")

tr <- TRUTH
u  <- as.numeric(q_true(0.90)); Su <- 1 - p_true(u)
sig_lim <- tr$sigma + tr$xi * (u - tr$mu); xi_lim <- tr$xi

# --- the truth's conditional law above u, in the exceedance variable z -------
# partial_moment_true integrates numerically, far too slow inside an optimiser.
# Two facts make it cheap. Beyond x = 8 the normal component is within 1e-12 of
# one, so the truth IS its GEV component there and gev_partial_moment is exact.
# Below that, tabulate once and interpolate with a monotone spline.
XCUT <- 8
.zt <- seq(0, XCUT - u, length.out = 800)
.pt <- vapply(.zt, function(zz) partial_moment_true(u + zz), numeric(1))
.sp <- splinefun(.zt, .pt, method = "hyman")
Fc    <- function(z) (p_true(u + z) - (1 - Su)) / Su
phic  <- function(z) {
  out <- numeric(length(z)); hi <- z > XCUT - u
  if (any(hi))  out[hi]  <- gev_partial_moment(u + z[hi], tr$mu, tr$sigma, tr$xi)
  if (any(!hi)) out[!hi] <- .sp(pmax(z[!hi], 0))
  out / Su
}
EZ    <- phic(0)
# E[(z-Z)^+] = z - E[Z] + E[(Z-z)^+], but the exceedances live on [0, Inf), so
# the lower partial moment is identically zero for z <= 0. Only the one-sided
# family evaluates it there (at t - c), and getting this wrong drives its
# estimating equation to a spurious boundary solution.
phiLc <- function(z) ifelse(z <= 0, 0, z - EZ + phic(pmax(z, 0)))

# --- identification functions g_p(t) under the truth -------------------------
g_quant <- function(t, p) p - Fc(t)
g_expec <- function(t, p) p * phic(t) - (1 - p) * phiLc(t)
g_elast <- function(s, a) function(t, p)
  (2 * a / s) * (p * phic(t) - (1 - p) * phiLc(t)) + (1 - a) * (p - Fc(t))
g_oneside <- function(cc) function(t, p)
  p * phic(t) - (1 - p) * (phiLc(t) - phiLc(t - cc))

# --- model functionals T(p | sigma, xi), threshold at zero -------------------
T_quant <- function(p, th) qgpd(p, 0, th[1], th[2])
T_expec <- function(p, th) gpd_expectile(p, 0, th[1], th[2])
T_elast <- function(s, a) function(p, th) gpd_elastile(p, 0, th[1], th[2], a, s)
T_oneside <- function(cc) function(p, th) gpd_onesided(p, 0, th[1], th[2], cc)

# --- solve (*) ---------------------------------------------------------------
solve_pt <- function(gfun, Tfun, w_fun, grid, start = c(sig_lim, xi_lim)) {
  p <- grid$p; wq <- grid$w_quad * w_fun(p)
  k <- wq > 0; p <- p[k]; wq <- wq[k]
  G <- function(th) {
    if (th[1] <= 0 || th[2] <= -0.4 || th[2] >= 0.9) return(c(1e6, 1e6))
    tv <- Tfun(p, th); if (any(!is.finite(tv))) return(c(1e6, 1e6))
    gv <- gfun(tv, p)
    h <- 1e-6 * pmax(abs(th), 1)             # dT/dtheta by central differences
    d1 <- (Tfun(p, th + c(h[1], 0)) - Tfun(p, th - c(h[1], 0))) / (2 * h[1])
    d2 <- (Tfun(p, th + c(0, h[2])) - Tfun(p, th - c(0, h[2]))) / (2 * h[2])
    c(sum(wq * gv * d1), sum(wq * gv * d2))
  }
  obj <- function(th) { v <- G(th); sum((v / c(1, 1))^2) }
  # several starts, keep the best root: the criterion has a flat ridge and a
  # single Nelder-Mead run can walk into the shape boundary instead
  starts <- list(start, c(1.2, 0.34), c(1.45, 0.21), c(1.35, 0.26), c(1.6, 0.18))
  best <- NULL
  for (st in starts) {
    o <- try(optim(st, obj, method = "Nelder-Mead",
                   control = list(reltol = 1e-14, maxit = 4000)), silent = TRUE)
    if (inherits(o, "try-error")) next
    o <- optim(o$par, obj, method = "Nelder-Mead",
               control = list(reltol = 1e-14, maxit = 4000))
    if (is.null(best) || o$value < best$value) best <- o
  }
  list(par = best$par, resid = sqrt(best$value))
}

grid <- make_level_grid(0, 24, 8)              # the fitting grid
W6 <- function(p) p^6

# elastile scale, fixed at the value the limiting GPD implies
s_el <- local({
  p <- grid$p; wq <- grid$w_quad * W6(p); k <- wq > 0; p <- p[k]; wq <- wq[k]
  te <- gpd_expectile(p, 0, sig_lim, xi_lim); tq <- qgpd(p, 0, sig_lim, xi_lim)
  L2 <- sum(wq * (p * phic(te) + (1 - p) * phiLc(te)))
  L1 <- sum(wq * (p * phic(tq) + (1 - p) * phiLc(tq)))
  L2 / L1
})
# the knot the fitter would use: k times the IQR of the exceedances themselves,
# i.e. of the truth's conditional law, not of the limiting GPD
Qc0 <- function(v) vapply(v, function(vv) as.numeric(q_true(1 - Su * (1 - vv))), numeric(1)) - u
cc  <- 4 * diff(Qc0(c(0.25, 0.75)))

EST <- list(
  "composite L1 (quantile)" = list(g_quant, T_quant),
  "composite L2 (expectile)" = list(g_expec, T_expec),
  "elastile a=0.5"          = list(g_elast(s_el, 0.5), T_elast(s_el, 0.5)),
  "one-sided k=4"           = list(g_oneside(cc), T_oneside(cc))
)

# --- exact references --------------------------------------------------------
kl <- function(par) {
  s <- exp(par[1]); xi <- par[2]
  f <- function(x) {
    d <- d_true(x)
    lg <- -log(s) - (1/xi + 1) * log1p(pmax(xi * (x - u)/s, -1 + 1e-12))
    o <- -d * lg / Su; o[!is.finite(o)] <- 0; o
  }
  v <- try(integrate(f, u, Inf, subdivisions = 4000, rel.tol = 1e-10)$value, silent = TRUE)
  if (inherits(v, "try-error") || !is.finite(v)) 1e10 else v
}
o <- optim(c(log(1.3), 0.3), kl, control = list(reltol = 1e-14, maxit = 3000))
mle_pt <- c(exp(o$par[1]), o$par[2])

# L-moments: lambda1 = int Q, lambda2 = int Q (2v-1), from the truth's conditional Q
gl <- gauss_legendre_01(2000)
Qc <- function(v) vapply(v, function(vv) as.numeric(q_true(1 - Su * (1 - vv))), numeric(1)) - u
Qv <- Qc(gl$x)
l1 <- sum(gl$w * Qv); l2 <- sum(gl$w * Qv * (2 * gl$x - 1))
xi_lm <- 2 - l1 / l2; s_lm <- l1 * (1 - xi_lm)
lmom_pt <- c(s_lm, xi_lm)

sink("out/pseudotrue-exact.txt", split = TRUE)
cat(sprintf("=== Pseudo-true GPD limits above u = %.4f, solved not simulated ===\n\n", u))
cat(sprintf("limiting GPD (threshold stability):   scale %.5f  shape %.5f\n", sig_lim, xi_lim))
cat(sprintf("elastile scale s = %.4f;  one-sided knot c = 4 x IQR = %.4f\n\n", s_el, cc))
cat("Conventional references, each solved exactly:\n")
cat(sprintf("  %-26s scale %.5f  shape %.5f   (KL minimiser)\n", "POT maximum likelihood",
            mle_pt[1], mle_pt[2]))
cat(sprintf("  %-26s scale %.5f  shape %.5f   (L-moment matching)\n\n", "POT L-moments",
            lmom_pt[1], lmom_pt[2]))
cat("Composite estimators, solving the estimating equation (*), w(p) = p^6:\n")
for (nm in names(EST)) {
  r <- solve_pt(EST[[nm]][[1]], EST[[nm]][[2]], W6, grid)
  cat(sprintf("  %-26s scale %.5f  shape %.5f   (|G| = %.1e)\n",
              nm, r$par[1], r$par[2], r$resid))
}
cat("\nWeight dependence, composite L2 (expectile), w(p) = p^m:\n")
for (m in c(0, 1, 2, 4, 6, 10, 20)) {
  r <- solve_pt(g_expec, T_expec, local({mm <- m; function(p) p^mm}), grid)
  cat(sprintf("  m = %-4d  scale %.5f  shape %.5f   (|G| = %.1e)\n",
              m, r$par[1], r$par[2], r$resid))
}
sink()
cat("\nwrote out/pseudotrue-exact.txt\n")

# ---------------------------------------------------------------------------
# The figure, redrawn from the exact limits. scripts/28 keeps the simulated
# version as an independent cross-check; these are the values to quote.
# ---------------------------------------------------------------------------
PT <- list(
  "POT MLE"        = mle_pt,
  "POT L-moments"  = lmom_pt,
  "composite L1"   = solve_pt(g_quant, T_quant, W6, grid)$par,
  "composite L2"   = solve_pt(g_expec, T_expec, W6, grid)$par,
  "one-sided k=4"  = solve_pt(g_oneside(cc), T_oneside(cc), W6, grid)$par
)
saveRDS(list(u = u, Su = Su, sig_lim = sig_lim, xi_lim = xi_lim, PT = PT),
        "out/pseudotrue-exact.rds")

xmax  <- u + qgpd(1 - 1e-3, 0, sig_lim, xi_lim)
xg    <- seq(u, xmax, length.out = 700)
Slim  <- sgpd(xg - u, 0, sig_lim, xi_lim)
Struth<- (1 - p_true(xg)) / Su
Sfit  <- lapply(PT, function(th) sgpd(xg - u, 0, th[1], th[2]))

set.seed(7)
ys <- r_true(N_OBS); ze <- sort(ys[ys > u] - u); m <- length(ze)
step_x <- u + c(0, ze[seq_len(m - 1)]); step_y <- c(1, (m - seq_len(m - 1)) / m)

nms  <- names(PT)
cols <- c("POT MLE" = "#a53a2b", "POT L-moments" = "#e67e22",
          "composite L1" = "#7d3c98", "composite L2" = "#197a45",
          "one-sided k=4" = "#1f5f8b")

png("out/fig-pseudotrue-tail.png", width = 1700, height = 700, res = 133)
par(mfrow = c(1, 2), mar = c(4.6, 4.8, 3.6, 1.2))
plot(NA, xlim = c(u, xmax), ylim = c(1e-3, 1), log = "y",
     xlab = "x", ylab = "conditional survival  P(X > x | X > u)",
     main = sprintf("Pseudo-true fits above u = %.2f", u))
lines(step_x, step_y, col = "grey62", lwd = 1.7, type = "s")
segments(u + ze[m - 1], 1/m, u + ze[m], 1/m, col = "grey62", lwd = 1.7)
points(u + ze[m], 1/m, col = "grey62", pch = 16, cex = 0.7)
lines(xg, Struth, col = "grey35", lwd = 4.4)
for (nm in nms) lines(xg, Sfit[[nm]], col = cols[nm], lwd = 2.3)
lines(xg, Slim, col = "black", lwd = 3, lty = 2)
legend("topright", c("truth (conditional)", sprintf("empirical, n = %d", N_OBS),
                     "limiting GPD", nms),
       col = c("grey35", "grey62", "black", cols[nms]),
       lwd = c(4.4, 1.7, 3, rep(2.3, length(nms))),
       lty = c(1, 1, 2, rep(1, length(nms))), bty = "n", cex = 0.72)

plot(NA, xlim = c(u, xmax), ylim = c(0.55, 2.4), log = "y",
     xlab = "x", ylab = "survival relative to the limiting GPD",
     main = "Distance from the limit")
abline(h = 1, col = "black", lwd = 3, lty = 2)
lines(xg, Struth / Slim, col = "grey35", lwd = 4.4)
for (nm in nms) lines(xg, Sfit[[nm]] / Slim, col = cols[nm], lwd = 2.3)
mtext("flat at 1 = converges to the limiting GPD", side = 1, line = 3.5,
      cex = 0.66, col = "grey35")
dev.off()
cat("wrote out/fig-pseudotrue-tail.png (exact limits)\n")

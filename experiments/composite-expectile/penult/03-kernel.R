# ---------------------------------------------------------------------------
# Effective threshold, and the kernel that produces it.
#
# EFFECTIVE THRESHOLD. An estimator fitted above u returns some shape xi*. The
# penultimate shape r'(x) is the shape a GPD would have if it were anchored at
# x. The level x_e where r'(x_e) = xi* is the level the estimator is behaving
# as if it had been handed -- its effective threshold -- reported as the return
# period R = 1/S(x_e). The reading is only honest if the whole GPD matches, not
# just the shape, so the penultimate scale at x_e carried back to u,
#     sigma_u = r(x_e) - r'(x_e) (x_e - u),
# is checked against the target's own scale.
#
# KERNEL. If xi* is a weighted average of local shapes, xi* = int K(s) r'(u+s) ds,
# then K says which part of the tail the estimator is reading. K is recovered
# without assuming the representation: perturb r by a smooth STEP at position s,
# which is a unit-mass bump in r' there, and K(s) = d xi*/d eps. A step in r is
# far better conditioned than a narrow bump in r', which would need a grid fine
# enough to resolve it.
#
# Two things fall out as checks. int K = 1 must hold for any Fisher-consistent
# estimator, since shifting r' by a constant c shifts the whole GPD's shape by c.
# And K >= 0 is a claim about the estimator, not a construction: an estimator
# that reads some part of the tail with the wrong sign would show it here.
# ---------------------------------------------------------------------------

## r'(x) = xi_star, first crossing above u, as a return period
effective_threshold <- function(base, xi_star, sigma_star = NA) {
  P <- base$P; u <- base$u; Su <- base$Su
  a <- base$a[-1]; yy <- base$y[-1]
  d <- P$rp(u + yy) - xi_star
  sgn <- which(d[-1] * d[-length(d)] < 0)
  if (!length(sgn)) return(list(R = NA_real_, x = NA_real_, sigma_u = NA_real_,
                                sigma_err = NA_real_, note = "no crossing above u"))
  i <- sgn[1]
  yof <- splinefun(base$a, base$y, method = "monoH.FC")
  f <- function(ee) P$rp(u + yof(ee)) - xi_star
  eo <- uniroot(f, c(a[i], a[i+1]), tol = 1e-13)$root
  xe <- u + yof(eo)
  sg_u <- P$r(xe) - P$rp(xe) * (xe - u)
  list(R = exp(eo) / Su, x = xe, sigma_u = sg_u,
       sigma_err = if (is.na(sigma_star)) NA_real_ else (sg_u - sigma_star)/sigma_star,
       note = if (length(sgn) > 1) sprintf("%d crossings", length(sgn)) else "")
}

## K(s) at a grid of exceedance levels e0, by perturbing r with a smooth step
kernel <- function(base, fit, e0 = seq(0.05, 13, length.out = 40), eps = 1e-3,
                   hwid = 0.4) {
  P <- base$P; u <- base$u
  ex0 <- excess_on(base)
  eb  <- ex0$e_of_y
  ss  <- function(v) { v <- pmin(1, pmax(0, v)); v*v*(3 - 2*v) }
  bt  <- fit(ex0)
  out <- vapply(e0, function(E0) {
    step <- function(x) ss((eb(x - u) - (E0 - hwid)) / (2 * hwid))
    a <- fit(excess_on(base, function(x) P$r(x) + eps * step(x)), start = bt)
    b <- fit(excess_on(base, function(x) P$r(x) - eps * step(x)), start = bt)
    (a[2] - b[2]) / (2 * eps)
  }, numeric(1))
  ## K is a density in s = y; its mass element on the level grid is K dy/de = K r
  yy <- splinefun(base$a, base$y, method = "monoH.FC")(e0)
  list(e0 = e0, y = yy, Su = base$Su, R = exp(e0) / base$Su, K = out,
       dens_e = out * P$r(u + yy), target = bt,
       total = sum(diff(e0) * (head(out * P$r(u + yy), -1) +
                               tail(out * P$r(u + yy), -1)) / 2))
}

## quantiles of K, as return periods. K is renormalised to integrate to one --
## int K = 1 is exact only if xi* really is a linear functional of r', so the
## measured total is reported alongside as how good that representation is.
kernel_quantiles <- function(kk, probs = c(.10, .25, .50, .75, .90)) {
  e <- kk$e0; m <- pmax(kk$dens_e, 0)
  cum <- c(0, cumsum(diff(e) * (head(m, -1) + tail(m, -1)) / 2)); cum <- cum / cum[length(cum)]
  vapply(probs, function(q) exp(approx(cum, e, q, ties = "ordered")$y) / kk$Su, numeric(1))
}

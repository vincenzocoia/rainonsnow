# ---------------------------------------------------------------------------
# Which weight functions are admissible, and how the two losses differ.
#
# The composite estimators minimise  R(theta) = int_0^1 w(p) r_p(theta) dp  with
#
#   r_p(theta) = E[ rho_p(Y, T(p | theta)) ]
#
# the population risk at level p. For the estimator to mean anything, R has to
# be finite, which constrains BOTH the weight's behaviour as p -> 1 and the
# tail of the truth. This script works both out, in closed form where possible
# and numerically as a check.
#
# Output: out/admissibility.txt, out/fig-admissibility.png
# ---------------------------------------------------------------------------
source("R/setup.R")
source("R/config.R")
dir.create("out", showWarnings = FALSE)
sink("out/admissibility.txt", split = TRUE)

# --- the two population risks, in closed form for a GEV truth ---------------
# quantile:  E[rho_p(Y,q)] = p E[(Y-q)^+] + (1-p) E[(q-Y)^+]
#                          = phi(q) + (1-p)(q - m)
risk_L1 <- function(p, q, mu, sigma, xi) {
  gev_partial_moment(q, mu, sigma, xi) + (1 - p) * (q - gev_mean(mu, sigma, xi))
}
# expectile: E[|p - I|(Y-e)^2] = p E[(Y-e)^2 1{Y>=e}] + (1-p) E[(e-Y)^2 1{Y<e}]
# and E[(e-Y)^2 1{Y<e}] = E[(Y-e)^2] - E[(Y-e)^2 1{Y>=e}].
risk_L2 <- function(p, e, mu, sigma, xi) {
  if (xi >= 0.5) return(rep(Inf, length(p)))
  m  <- gev_mean(mu, sigma, xi)
  v  <- sigma^2 * (gamma(1 - 2 * xi) - gamma(1 - xi)^2) / xi^2   # variance
  up <- gev_second_partial_moment(e, mu, sigma, xi)
  tot <- v + (m - e)^2                                           # E[(Y-e)^2]
  p * up + (1 - p) * (tot - up)
}

cat("=== 1. Closed form for psi(x) = E[(X-x)^2 1{X>=x}] ===\n")
cat("(a) at the lower endpoint psi must equal E[(X - lo)^2] = Var + (mean - lo)^2;\n")
cat("(b) psi'(x) = -2 phi(x) identically, checked by central difference against\n")
cat("    the independently derived first partial moment;\n")
cat("(c) Monte Carlo, which is only a weak check for large xi -- E[(X-x)^4] is\n")
cat("    infinite once xi >= 1/4, so the MC estimator itself has infinite variance.\n\n")
for (par in list(c(0, 1, 0.2), c(0, 1, 0.35), c(0, 1, 0.45), c(2, 3, 0.1), c(0, 1, -0.2))) {
  mu <- par[1]; sigma <- par[2]; xi <- par[3]
  lo <- gev_endpoints(mu, sigma, xi)[1]
  m  <- gev_mean(mu, sigma, xi)
  v  <- sigma^2 * (gamma(1 - 2 * xi) - gamma(1 - xi)^2) / xi^2
  cat(sprintf("GEV(%.0f,%.0f,%+.2f)\n", mu, sigma, xi))
  if (is.finite(lo))
    cat(sprintf("   (a) psi(lo) = %.8f    Var + (m-lo)^2 = %.8f\n",
                gev_second_partial_moment(lo, mu, sigma, xi), v + (m - lo)^2))
  xs <- qgev(c(0.3, 0.5, 0.9, 0.99), mu, sigma, xi)
  h  <- 1e-5 * pmax(1, abs(xs))
  d  <- (gev_second_partial_moment(xs + h, mu, sigma, xi) -
         gev_second_partial_moment(xs - h, mu, sigma, xi)) / (2 * h)
  cat(sprintf("   (b) max rel error in psi'(x) = -2 phi(x): %.2e\n",
              max(abs(d + 2 * gev_partial_moment(xs, mu, sigma, xi)) /
                  abs(2 * gev_partial_moment(xs, mu, sigma, xi)))))
  set.seed(4); y <- qgev(runif(4e6), mu, sigma, xi)
  mc <- sapply(xs, function(x) mean((y - x)^2 * (y >= x)))
  cat(sprintf("   (c) MC / closed form at Q(.3, .5, .9, .99): %s\n",
      paste(sprintf("%.3f", mc / gev_second_partial_moment(xs, mu, sigma, xi)), collapse = "  ")))
}

cat("\n=== 2. How each risk decays as p -> 1 ===\n")
cat("Write u = 1 - p. With the model at the truth, the tail behaviour is\n\n")
cat("   quantile :  r_p ~ u^(1 - xi)      (from the (1-p)(q - m) term)\n")
cat("   expectile:  r_p ~ u^(1 - 2 xi)    (from the (1-p) E[(e-Y)^2] term)\n\n")
cat("so int^1 w(p) r_p dp converges for any BOUNDED w whenever xi < 1 (quantile)\n")
cat("or xi < 1/2 (expectile). Measured exponents (slope of log r vs log u):\n\n")
for (xi in c(0.1, 0.2, 0.3, 0.45)) {
  u <- 10^seq(-6, -3, length.out = 40); p <- 1 - u
  q <- qgev(p, 0, 1, xi); e <- gev_expectile(p, 0, 1, xi)
  s1 <- coef(lm(log(risk_L1(p, q, 0, 1, xi)) ~ log(u)))[2]
  s2 <- coef(lm(log(risk_L2(p, e, 0, 1, xi)) ~ log(u)))[2]
  cat(sprintf("xi = %.2f   quantile slope %6.3f (predicted %5.2f)   expectile slope %6.3f (predicted %5.2f)\n",
              xi, s1, 1 - xi, s2, 1 - 2 * xi))
}

cat("\n=== 3. The condition on an UNBOUNDED weight ===\n")
cat("For w(p) ~ (1 - p)^(-a) as p -> 1, int w(p) r_p dp converges iff\n\n")
cat("   quantile :  a < 2 - xi\n")
cat("   expectile:  a < 2 - 2 xi\n\n")
cat("At xi = 0.2 that is a < 1.8 and a < 1.6. So a weight that blows up like\n")
cat("1/(1-p) (a = 1) is still admissible for both, and one like 1/(1-p)^2 is\n")
cat("admissible for neither. The weight used here is bounded -- it is exactly 1\n")
cat("on [p1, 1) and exactly 0 below p0 -- so it satisfies both conditions with\n")
cat("room to spare, for every xi in the range the estimators explore.\n")

cat("\n=== 4. The condition that actually separates L1 from L2 ===\n")
cat("The quantile risk needs E|Y| < Inf (xi < 1). The expectile risk needs\n")
cat("E[Y^2] < Inf (xi < 1/2): above that, E[(e - Y)^2 1{Y < e}] is infinite for\n")
cat("every e and the whole criterion is +Inf regardless of the weight.\n\n")
for (xi in c(0.3, 0.45, 0.5, 0.6)) {
  q <- qgev(0.99, 0, 1, xi)
  r1 <- risk_L1(0.99, q, 0, 1, xi)
  e <- gev_expectile(0.99, 0, 1, xi)
  r2 <- risk_L2(0.99, e, 0, 1, xi)
  cat(sprintf("xi = %.2f   r_0.99 quantile = %8.4f   expectile = %s\n",
              xi, r1, if (is.finite(r2)) sprintf("%.4f", r2) else "Inf"))
}
cat("\nThis narrows where the L2 version can be used at all, relative to L1, and\n")
cat("it narrows it exactly in the heavy-tailed direction the method is aimed at.\n")
cat("The fix is standard: minimise the loss DIFFERENCE against a fixed reference,\n")
cat("sum_i int w(p) [rho_p(y_i, T(p)) - rho_p(y_i, T_ref(p))] dp, which has the\n")
cat("same minimiser and is finite whenever E|Y| < Inf. In this experiment\n")
cat("xi = 0.2, so the raw criterion is finite and nothing was needed.\n")

cat("\n=== 5. Finiteness of the criterion actually used here ===\n")
tail_share <- function(p_from) {
  k <- GRID$p > p_from
  sum(GRID$w_quad[k] * W_FUN(GRID$p[k])) / sum(GRID$w_quad * W_FUN(GRID$p))
}
cat(sprintf("weight support: [%g, 1), bounded by 1\n", P0))
cat(sprintf("int w(p) dp = %.6f\n", sum(GRID$w_quad * W_FUN(GRID$p))))
for (pf in c(0.99, 0.999, 0.9999, 0.99999))
  cat(sprintf("share of the weight above p = %-8g : %.4f\n", pf, tail_share(pf)))
xi <- TRUTH$xi
pg <- GRID$p
cat(sprintf("\npopulation risk of the weighted criterion at the true GEV tail (xi = %g):\n", xi))
cat(sprintf("  quantile  = %.6f\n", sum(GRID$w_quad * W_FUN(pg) *
      risk_L1(pg, qgev(pg, 0, 1, xi), 0, 1, xi))))
cat(sprintf("  expectile = %.6f\n", sum(GRID$w_quad * W_FUN(pg) *
      risk_L2(pg, gev_expectile(pg, 0, 1, xi), 0, 1, xi))))
cat("Both finite, as the conditions above require.\n")
sink()

png("out/fig-admissibility.png", width = 1500, height = 620, res = 133)
par(mfrow = c(1, 2), mar = c(4.4, 4.6, 3.2, 1.2))
u <- 10^seq(-7, -1.4, length.out = 200); p <- 1 - u
cols <- c("#2980b9", "#1e8449", "#c0392b", "#8e44ad")
xis <- c(0.1, 0.2, 0.3, 0.45)
for (k in 1:2) {
  r <- sapply(seq_along(xis), function(i) { xi <- xis[i]
    if (k == 1) risk_L1(p, qgev(p, 0, 1, xi), 0, 1, xi)
    else        risk_L2(p, gev_expectile(p, 0, 1, xi), 0, 1, xi) })
  plot(u, r[, 1], type = "n", log = "xy", ylim = range(r),
       xlab = "1 - p", ylab = "population risk at level p",
       main = sprintf("%s risk as p -> 1", c("Quantile (L1)", "Expectile (L2)")[k]))
  for (i in seq_along(xis)) lines(u, r[, i], col = cols[i], lwd = 2.4)
  legend("bottomright", sprintf("xi = %.2f  (slope %.2f)", xis,
         if (k == 1) 1 - xis else 1 - 2 * xis), col = cols, lwd = 2.4, bty = "n", cex = 0.8)
}
dev.off()
cat("\nwrote out/admissibility.txt and out/fig-admissibility.png\n")

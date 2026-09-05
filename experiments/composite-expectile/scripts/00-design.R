# ---------------------------------------------------------------------------
# Design of the experiment.
#
# Fixes the data-generating process and the weight function, and works out what
# each estimator converges to (its "pseudo-true" parameter) by fitting it to one
# very large sample. That asymptotic target is the bias each estimator carries
# no matter how much data it is given.
#
# Outputs: out/design.txt, out/fig-design.png
# ---------------------------------------------------------------------------
source("R/setup.R")
source("R/config.R")
dir.create("out", showWarnings = FALSE)
sink("out/design.txt", split = TRUE)

tr <- TRUTH
cat("=== Data-generating process ===\n")
cat("F(x) = F_GEV(x; mu, sigma, xi) * Phi((x - a)/b)   [ = cdf of max(GEV, Normal) ]\n")
cat(sprintf("  GEV:    mu = %g, sigma = %g, xi = %g\n", tr$mu, tr$sigma, tr$xi))
cat(sprintf("  Normal: a  = %g, b     = %g\n", tr$a, tr$b))
cat(sprintf("  true support: [%g, Inf)\n", gev_endpoints(tr$mu, tr$sigma, tr$xi)[1]))
cat("The normal contaminates the body and vanishes in the tail, so the GEV is\n")
cat("the exactly correct TAIL model and a wrong BODY model.\n\n")

cat("=== Where the contamination dies out ===\n")
cat("contamination(p) = P(Normal > Q_true(p)): how much pull the normal still has.\n\n")
pp <- c(0.5, 0.8, 0.9, 0.95, 0.97, 0.98, 0.99, 0.995, 0.999)
qt <- q_true(pp); qg <- qgev(pp, tr$mu, tr$sigma, tr$xi)
tab <- data.frame(p = pp, return_period = 1 / (1 - pp), Q_true = qt, Q_GEV = qg,
                  pct_excess = 100 * (qt - qg) / qg,
                  contamination = contamination(qt), w = W_FUN(pp))
print(format(tab, digits = 4), row.names = FALSE)

cat(sprintf("\n=== Weight function ===\nsmoothstep, 0 below p0 = %g, 1 above p1 = %g\n", P0, P1))
cat(sprintf("p0 sits where contamination ~ 1e-2 (quantiles still %.1f%% above the GEV),\n",
            100 * (q_true(P0) - qgev(P0, tr$mu, tr$sigma, tr$xi)) / qgev(P0, tr$mu, tr$sigma, tr$xi)))
cat(sprintf("p1 where contamination ~ 1e-4 (quantiles %.2f%% above the GEV).\n",
            100 * (q_true(P1) - qgev(P1, tr$mu, tr$sigma, tr$xi)) / qgev(P1, tr$mu, tr$sigma, tr$xi)))
cat(sprintf("Quadrature: %d levels (%d panels x %d Gauss-Legendre nodes) on [%g, 1).\n",
            GRID$n_nodes, N_PANEL, N_GL, P0))

cat("\n=== Asymptotic targets (fit to one sample of n = 2e6) ===\n")
cat("What each estimator converges to. Any gap from the true GEV (0, 1, 0.2) is\n")
cat("bias that no sample size removes.\n\n")
set.seed(20260905)
big <- r_true(2e6)
pseudo <- rbind(
  MLE      = fit_mle(big),
  Lmoments = fit_lmom(big),
  CQE_L1   = fit_cqe(big, GRID, W_FUN),
  CEE_L2   = fit_cee(big, GRID, W_FUN)
)
colnames(pseudo) <- c("mu", "sigma", "xi")
print(round(pseudo, 4))
cat(sprintf("\ntrue GEV tail parameters: mu = %g, sigma = %g, xi = %g\n",
            tr$mu, tr$sigma, tr$xi))

cat("\n--- implied return levels at the asymptotic targets ---\n")
Tp <- c(2, 5, 10, 20, 50, 100, 200, 500, 1000)
rl <- rbind(truth = q_true(1 - 1 / Tp),
            t(apply(pseudo, 1, function(th) qgev(1 - 1 / Tp, th[1], th[2], th[3]))))
colnames(rl) <- paste0("T=", Tp)
print(round(rl, 3))
cat("\n--- asymptotic % bias in return level ---\n")
print(round(100 * (rl[-1, , drop = FALSE] -
                     rep(rl[1, ], each = nrow(rl) - 1)) / rep(rl[1, ], each = nrow(rl) - 1), 2))
sink()

# --- design figure ---------------------------------------------------------
png("out/fig-design.png", width = 1500, height = 520, res = 130)
par(mfrow = c(1, 3), mar = c(4.2, 4.2, 3, 1), cex.main = 1.0)
xs <- seq(-3, 14, length.out = 800)
plot(xs, p_true(xs), type = "l", lwd = 2.5, col = "black", ylim = c(0, 1),
     xlab = "x", ylab = "cdf", main = "Truth vs its GEV tail")
lines(xs, pgev(xs, tr$mu, tr$sigma, tr$xi), lwd = 2, col = "#c0392b", lty = 2)
lines(xs, pnorm(xs, tr$a, tr$b), lwd = 1.5, col = "#2980b9", lty = 3)
legend("bottomright", c("truth (product)", "GEV component", "normal component"),
       col = c("black", "#c0392b", "#2980b9"), lty = c(1, 2, 3), lwd = 2, bty = "n", cex = 0.8)

pg <- seq(0.5, 0.9999, length.out = 600)
qtg <- q_true(pg); qgg <- qgev(pg, tr$mu, tr$sigma, tr$xi)
plot(1 / (1 - pg), 100 * (qtg - qgg) / qgg, type = "l", lwd = 2.5, log = "x",
     xlab = "return period (years)", ylab = "% by which truth exceeds the GEV",
     main = "Model error left by the\ncontamination")
abline(h = 0, col = "grey70"); abline(v = 1 / (1 - c(P0, P1)), col = "#27ae60", lty = c(2, 3))
text(1 / (1 - P0), par("usr")[4] * 0.75, " p0", col = "#27ae60", adj = 0, cex = 0.9)
text(1 / (1 - P1), par("usr")[4] * 0.55, " p1", col = "#27ae60", adj = 0, cex = 0.9)

plot(1 / (1 - pg), W_FUN(pg), type = "l", lwd = 2.5, col = "#27ae60", log = "x",
     ylim = c(0, 1.05), xlab = "return period (years)", ylab = "weight",
     main = "Weight function w(p)")
lines(1 / (1 - pg), contamination(qtg) / max(contamination(qtg)), lwd = 2,
      col = "#8e44ad", lty = 2)
legend("right", c("w(p)", "contamination (scaled)"), col = c("#27ae60", "#8e44ad"),
       lty = c(1, 2), lwd = 2, bty = "n", cex = 0.8)
dev.off()
cat("\nwrote out/design.txt and out/fig-design.png\n")

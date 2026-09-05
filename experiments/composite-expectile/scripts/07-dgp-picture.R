# ---------------------------------------------------------------------------
# What the data-generating process actually looks like.
#
# Densities of the two components and of their maximum, and the
# frequency-magnitude curve (return period vs return level) of all three. The
# frequency-magnitude panel is the one that shows whether the normal is really
# contaminating anything: where the truth sits above its GEV component, it is.
#
# Two settings are drawn: the one used in the experiment, and the heavier
# alternative Vincenzo suggested, GEV(23, 22.2, 0.39) with N(61.9, 16.9).
#
# Output: out/fig-dgp.png, out/dgp-comparison.txt
# ---------------------------------------------------------------------------
source("R/setup.R")
source("R/config.R")
dir.create("out", showWarnings = FALSE)

SETTINGS <- list(
  `experiment: GEV(0, 1, 0.20) v N(1.5, 0.8)` =
    list(mu = 0, sigma = 1, xi = 0.20, a = 1.5, b = 0.8),
  `suggested: GEV(23, 22.2, 0.39) v N(61.9, 16.9)` =
    list(mu = 23, sigma = 22.2, xi = 0.39, a = 61.9, b = 16.9)
)

Tgrid <- exp(seq(log(1.05), log(2000), length.out = 400))
pgrid <- 1 - 1 / Tgrid

png("out/fig-dgp.png", width = 1600, height = 1080, res = 132)
par(mfrow = c(2, 2), mar = c(4.3, 4.4, 3.2, 1.2))
for (nm in names(SETTINGS)) {
  tr <- SETTINGS[[nm]]
  lo <- gev_endpoints(tr$mu, tr$sigma, tr$xi)[1]
  hi <- q_true(0.995, tr)
  xs <- seq(max(lo, tr$a - 4 * tr$b) - 0.02 * (hi - lo), hi, length.out = 1200)
  xs <- xs[xs > lo]

  # --- densities ---
  dG <- dgev(xs, tr$mu, tr$sigma, tr$xi)
  dN <- dnorm(xs, tr$a, tr$b)
  dT <- d_true(xs, tr)
  plot(xs, dT, type = "l", lwd = 2.8, col = "black", ylim = c(0, max(dT, dG, dN)),
       xlab = "x", ylab = "density", main = paste0("Densities\n", nm), cex.main = 0.9)
  lines(xs, dG, lwd = 2, col = "#c0392b", lty = 2)
  lines(xs, dN, lwd = 2, col = "#2980b9", lty = 3)
  legend("topright", c("truth = max(GEV, N)", "GEV component", "normal component"),
         col = c("black", "#c0392b", "#2980b9"), lty = c(1, 2, 3), lwd = 2,
         bty = "n", cex = 0.75)

  # --- frequency-magnitude ---
  rlT <- q_true(pgrid, tr)
  rlG <- qgev(pgrid, tr$mu, tr$sigma, tr$xi)
  rlN <- qnorm(pgrid, tr$a, tr$b)
  plot(Tgrid, rlT, type = "l", lwd = 2.8, col = "black", log = "x",
       ylim = range(c(rlT, rlN[pgrid > 0.05]), finite = TRUE),
       xlab = "return period T (years)", ylab = "return level",
       main = paste0("Frequency-magnitude\n", nm), cex.main = 0.9)
  lines(Tgrid, rlG, lwd = 2, col = "#c0392b", lty = 2)
  lines(Tgrid, rlN, lwd = 2, col = "#2980b9", lty = 3)
  legend("topleft", c("truth = max(GEV, N)", "GEV component", "normal component"),
         col = c("black", "#c0392b", "#2980b9"), lty = c(1, 2, 3), lwd = 2,
         bty = "n", cex = 0.75)
}
dev.off()

sink("out/dgp-comparison.txt", split = TRUE)
for (nm in names(SETTINGS)) {
  tr <- SETTINGS[[nm]]
  cat("\n=== ", nm, " ===\n", sep = "")
  pp <- c(0.5, 0.8, 0.9, 0.95, 0.98, 0.99, 0.995, 0.999)
  qt <- q_true(pp, tr); qg <- qgev(pp, tr$mu, tr$sigma, tr$xi)
  qn <- qnorm(pp, tr$a, tr$b)
  et <- expectile_true(pp, tr); eg <- gev_expectile(pp, tr$mu, tr$sigma, tr$xi)
  print(format(data.frame(
    p = pp, `T` = 1 / (1 - pp),
    Q_truth = qt, Q_GEV = qg, Q_normal = qn,
    `Q_pct_excess` = 100 * (qt - qg) / qg,
    `E_pct_excess` = 100 * (et - eg) / eg,
    `P(N>Q_truth)` = contamination(qt, tr),
    check.names = FALSE), digits = 4), row.names = FALSE)
  cat(sprintf("\nwhich component wins at the median: GEV Q(.5) = %.2f, normal Q(.5) = %.2f\n",
              qgev(0.5, tr$mu, tr$sigma, tr$xi), qnorm(0.5, tr$a, tr$b)))
  cat(sprintf("mean of truth = %.3f, mean of GEV component = %.3f (ratio %.2f)\n",
              mean_true(tr), gev_mean(tr$mu, tr$sigma, tr$xi),
              mean_true(tr) / gev_mean(tr$mu, tr$sigma, tr$xi)))
}
sink()
cat("\nwrote out/fig-dgp.png and out/dgp-comparison.txt\n")

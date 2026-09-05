# ---------------------------------------------------------------------------
# The same comparison on the heavier setting Vincenzo suggested:
# GEV(23, 22.2, 0.39) contaminated by N(61.9, 16.9).
#
# Two things make it a useful stress test. The shape is 0.39, not far below the
# xi = 1/2 boundary where the raw expectile criterion stops being finite. And the
# contamination is milder on the quantile scale (2.5% at T = 10, against 11% in
# the main setting), so the tail-focused estimators have less bias to win back.
#
# Output: out/fits-altdgp.rds, out/alt-dgp.txt, out/fig-alt-dgp.png
# ---------------------------------------------------------------------------
source("R/setup.R")
source("R/config.R")
library(parallel)

TRUTH <- list(mu = 23, sigma = 22.2, xi = 0.39, a = 61.9, b = 16.9)
P0_SET <- c(0.5, 0.95)   # trimmed: the mixture study in scripts/11 is the better use of time
n <- N_OBS
set.seed(9090)
datasets <- lapply(seq_len(N_REP), function(i) r_true(n, TRUTH))

wfun <- function(p0) if (p0 <= 0) function(p) rep(1, length(p))
                     else make_weight(p0, min(p0 + 0.03, 0.999))

cat("baselines ... "); t0 <- Sys.time()
base <- mclapply(datasets, function(y) { l <- fit_lmom(y)
  list(lmom = l, mle = fit_mle(y, start = l)) }, mc.cores = detectCores())
fits <- list(mle = do.call(rbind, lapply(base, `[[`, "mle")),
             lmom = do.call(rbind, lapply(base, `[[`, "lmom")))
cat(sprintf("%.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))

for (p0 in P0_SET) {
  g <- make_level_grid(p0, N_PANEL, N_GL); wf <- wfun(p0)
  for (ty in c("quantile", "expectile")) {
    cat(sprintf("p0 = %.2f  %-9s ... ", p0, ty)); t0 <- Sys.time()
    fits[[sprintf("%s|%.2f", ty, p0)]] <- do.call(rbind, mclapply(datasets,
      function(y) fit_composite(y, g, wf, ty, start = fit_lmom(y)),
      mc.cores = detectCores()))
    cat(sprintf("%.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  }
}
saveRDS(list(truth = TRUTH, fits = fits, n = n), "out/fits-altdgp.rds")

Tp <- RETURN_PERIODS; truth_rl <- q_true(1 - 1 / Tp, TRUTH)
summ <- function(m) {
  rl <- t(apply(m, 1, function(th) if (any(is.na(th))) rep(NA_real_, length(Tp))
                else qgev(1 - 1 / Tp, th[1], th[2], th[3])))
  rl <- rl[complete.cases(rl), , drop = FALSE]
  e <- sweep(rl, 2, truth_rl, "-")
  list(bias = colMeans(e), sd = apply(rl, 2, sd), mse = colMeans(e^2),
       xi = median(m[, 3], na.rm = TRUE))
}
S <- lapply(fits, summ)
Tshow <- c(10, 20, 50, 100, 200, 500, 1000)
idx <- sapply(Tshow, function(t) which.min(abs(Tp - t)))

sink("out/alt-dgp.txt", split = TRUE)
cat(sprintf("=== Heavier setting: GEV(%g, %g, %g) contaminated by N(%g, %g) ===\n",
            TRUTH$mu, TRUTH$sigma, TRUTH$xi, TRUTH$a, TRUTH$b))
cat(sprintf("n = %d, %d replicates\n\n", n, N_REP))
cat("--- median fitted shape (true xi = 0.39) ---\n")
print(round(sapply(S, `[[`, "xi"), 3))
for (lab in c("mse ratio vs GEV MLE", "bias", "sd")) {
  cat(sprintf("\n--- %s ---\n", lab))
  tab <- t(sapply(S, function(s) switch(lab, `mse ratio vs GEV MLE` = s$mse / S$mle$mse,
                                        bias = s$bias, sd = s$sd)[idx]))
  colnames(tab) <- paste0("T=", round(Tp[idx])); print(round(tab, 3))
}
sink()

png("out/fig-alt-dgp.png", width = 1600, height = 660, res = 134)
par(mfrow = c(1, 2), mar = c(4.4, 4.6, 3.2, 1.2))
cols <- colorRampPalette(c("#d5d8dc", "#154360"))(length(P0_SET))
for (ty in c("quantile", "expectile")) {
  plot(Tp, rep(1, length(Tp)), type = "n", log = "xy", ylim = c(0.3, 60),
       xlab = "return period T (years)", ylab = "MSE relative to GEV MLE",
       main = sprintf("Composite %s (%s)\nGEV(23, 22.2, 0.39) v N(61.9, 16.9)", ty,
                      if (ty == "quantile") "L1" else "L2"), cex.main = 0.9)
  abline(h = 1, col = "grey50", lwd = 1.5)
  lines(Tp, S$lmom$mse / S$mle$mse, col = "#e67e22", lwd = 3, lty = 2)
  for (i in seq_along(P0_SET))
    lines(Tp, S[[sprintf("%s|%.2f", ty, P0_SET[i])]]$mse / S$mle$mse, col = cols[i], lwd = 2.6)
  legend("topright", c(sprintf("p0 = %.2f", P0_SET), "GEV L-moments"),
         col = c(cols, "#e67e22"), lwd = 2.6,
         lty = c(rep(1, length(P0_SET)), 2), bty = "n", cex = 0.7)
}
dev.off()
cat("\nwrote out/fig-alt-dgp.png and out/alt-dgp.txt\n")

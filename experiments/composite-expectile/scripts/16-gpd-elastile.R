# ---------------------------------------------------------------------------
# The alpha-elastile inside the GPD study.
#
# Same setup as scripts/15: a three-parameter GPD fitted to the whole sample with
# a weight rising from p = 0, smooth grafted onto an empirical body with that
# same weight. Here the loss is swept from pure quantile (alpha = 0) to pure
# expectile (alpha = 1). The reference stays POT-MLE at the 90th percentile,
# hard grafted -- the practitioner's method.
#
# Uses R/graft_fast.R, which evaluates the graft directly instead of asking
# distplyr to recompute the correction integral at every point; it agrees with
# distplyr to 1e-10 when both are given the same analytic weight derivative and
# is about 200 times faster.
#
# Output: out/fits-gpd-elastile.rds, out/gpd-elastile.txt, out/fig-gpd-elastile.png
# ---------------------------------------------------------------------------
.libPaths(c("/home/user/Rlib-graft", .libPaths()))
source("R/setup.R"); source("R/config.R")
source("R/gpd.R"); source("R/gpd_estimators.R")
source("R/smoothgraft.R"); source("R/graft_fast.R")
suppressMessages({library(distionary); library(distplyr)})
library(parallel)

n <- N_OBS; NC <- detectCores()
Tp <- RETURN_PERIODS; ex <- 1 / Tp; truth_rl <- q_true(1 - ex)
set.seed(4242 + n)
datasets <- lapply(seq_len(N_REP), function(i) r_true(n))

W  <- function(p) p^6                       # the best weight from scripts/15
WD <- function(p) 6 * p^5
grid0 <- make_level_grid(0, N_PANEL, N_GL)
ALPHAS <- c(0, 0.1, 0.25, 0.5, 0.75, 0.9, 1)

RL <- list()
RL[["POT-MLE u=0.90"]] <- do.call(rbind, mclapply(datasets,
  function(y) pot_return_levels(y, fit_pot_mle(y, 0.90), ex), mc.cores = NC))
RL[["POT-Lmom u=0.90"]] <- do.call(rbind, mclapply(datasets,
  function(y) pot_return_levels(y, fit_pot_lmom(y, 0.90), ex), mc.cores = NC))
PAR <- list()
for (a in ALPHAS) {
  nm <- sprintf("alpha=%.2f", a)
  cat(sprintf("%s fit ... ", nm)); t0 <- Sys.time()
  PAR[[nm]] <- do.call(rbind, mclapply(datasets,
    function(y) fit_gpd_elastile(y, grid0, W, a), mc.cores = NC))
  cat(sprintf("%.1f min, graft ... ", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  t0 <- Sys.time()
  RL[[paste0(nm, " + graft")]] <- do.call(rbind, mclapply(seq_along(datasets),
    function(i) as.numeric(graft_fast_return_levels(datasets[[i]], PAR[[nm]][i, ],
                                                    W, ex, "gpd", w_deriv = WD)),
    mc.cores = NC))
  cat(sprintf("%.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
saveRDS(list(RL = RL, PAR = PAR, alphas = ALPHAS, Tp = Tp, truth = truth_rl),
        "out/fits-gpd-elastile.rds")

REF <- RL[["POT-MLE u=0.90"]]
summ <- function(r) {
  ok <- complete.cases(r) & complete.cases(REF)
  e <- sweep(r[ok, , drop = FALSE], 2, truth_rl, "-")
  er <- sweep(REF[ok, , drop = FALSE], 2, truth_rl, "-")
  list(ratio = colMeans(e^2) / colMeans(er^2),
       se = apply(e^2 - er^2, 2, sd) / sqrt(nrow(e)) / colMeans(er^2),
       medae = apply(abs(e), 2, median), win = colMeans(abs(e) < abs(er)),
       bias = colMeans(e), sd = apply(r[ok, , drop = FALSE], 2, sd), nfail = sum(!ok))
}
S <- lapply(RL, summ)
Tshow <- c(2, 10, 20, 50, 100, 200, 500, 1000)
idx <- sapply(Tshow, function(t) which.min(abs(Tp - t)))

sink("out/gpd-elastile.txt", split = TRUE)
cat(sprintf("=== alpha-elastile GPD + smooth graft, n = %d, %d replicates ===\n", n, N_REP))
cat("Weight w(p) = p^6 for both fitting and handover. All ratios against\n")
cat("POT-MLE(0.90) hard grafted onto an empirical body.\n\n")
cat("--- median fitted shape (GPD tail xi ~ 0.2) ---\n")
print(round(sapply(PAR, function(m) median(m[, 3], na.rm = TRUE)), 3))
for (lab in c("ratio", "se", "medae", "win", "bias", "sd")) {
  cat(sprintf("\n--- %s ---\n", switch(lab, ratio = "MSE relative to POT-MLE(0.90)",
    se = "MC standard error of that ratio", medae = "median absolute error",
    win = "fraction closer to the truth than POT-MLE(0.90)", bias = "bias", sd = "sd")))
  tab <- t(sapply(S, function(s) s[[lab]][idx])); colnames(tab) <- paste0("T=", round(Tp[idx]))
  print(round(tab, 3))
}
cat("\n--- alpha with the lowest MSE at each return period ---\n")
Sa <- S[paste0(sprintf("alpha=%.2f", ALPHAS), " + graft")]
best <- sapply(idx, function(k) ALPHAS[which.min(sapply(Sa, function(s) s$ratio[k]))])
print(setNames(best, paste0("T=", round(Tp[idx]))))
cat("--- and relative to the better of alpha = 0 and alpha = 1 ---\n")
gain <- sapply(idx, function(k) { m <- sapply(Sa, function(s) s$ratio[k])
  min(m) / min(m[1], m[length(m)]) })
print(round(setNames(gain, paste0("T=", round(Tp[idx]))), 3))
cat("\nfailed variants:\n"); print(sapply(S, `[[`, "nfail"))
sink()

png("out/fig-gpd-elastile.png", width = 1650, height = 640, res = 133)
par(mfrow = c(1, 3), mar = c(4.4, 4.6, 3.4, 1.2))
cols <- colorRampPalette(c("#2980b9", "#7d3c98", "#1e8449"))(length(ALPHAS))
plot(Tp, rep(1, length(Tp)), type = "n", log = "xy", ylim = c(0.1, 5),
     xlab = "return period T", ylab = "MSE relative to POT-MLE(0.90)",
     main = sprintf("GPD elastile + graft (n = %d)", n))
abline(h = 1, col = "grey50", lwd = 1.5)
lines(Tp, S[["POT-Lmom u=0.90"]]$ratio, col = "#e67e22", lwd = 3, lty = 2)
for (i in seq_along(ALPHAS)) lines(Tp, Sa[[i]]$ratio, col = cols[i], lwd = 2.4)
legend("bottomleft", c(sprintf("alpha = %.2f", ALPHAS), "POT-Lmom"),
       col = c(cols, "#e67e22"), lwd = 2.4, lty = c(rep(1, length(ALPHAS)), 2),
       bty = "n", cex = 0.65)
for (k in c(which.min(abs(Tp - 100)), which.min(abs(Tp - 1000)))) {
  r <- sapply(Sa, function(s) s$ratio[k])
  plot(ALPHAS, r, type = "b", pch = 19, lwd = 2.4, col = "#7d3c98",
       ylim = range(c(r, 1, S[["POT-Lmom u=0.90"]]$ratio[k])),
       xlab = "alpha  (0 = quantile, 1 = expectile)",
       ylab = "MSE relative to POT-MLE(0.90)", main = sprintf("T = %.0f", Tp[k]))
  abline(h = 1, col = "#c0392b", lwd = 2)
  abline(h = S[["POT-Lmom u=0.90"]]$ratio[k], col = "#e67e22", lwd = 2, lty = 2)
  legend("topright", c("elastile + graft", "POT-MLE(0.90)", "POT-Lmom(0.90)"),
         col = c("#7d3c98", "#c0392b", "#e67e22"), lwd = 2.4, lty = c(1,1,2),
         bty = "n", cex = 0.7)
}
dev.off()
cat("\nwrote out/fig-gpd-elastile.png and out/gpd-elastile.txt\n")

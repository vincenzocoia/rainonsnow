# ---------------------------------------------------------------------------
# 19. The one-sided inverted-Huber M-quantile inside the GPD study.
#
# psi(u) = max(u, -c): fully linear above the fitted level, capped at -c below
# it. Section 18 shows this removes the mean anchoring exactly -- at matched
# effective level the surviving body contamination is 0.000% against the
# quantile's 0.569% and the expectile's 15.7%, while the upside influence stays
# unbounded. The question here is the one that sank the plug-in fix in section
# 10: does it survive n = 100?
#
# Same protocol as scripts/16 so the numbers are directly comparable: three-
# parameter GPD, weight p^6 for both fitting and handover, smooth grafted onto
# an empirical body, ratios against POT-MLE(0.90) hard grafted. c = k * IQR(y),
# fixed from the data before optimising, so rho does not move with theta.
#
# The symmetric inversion psi(u) = sign(u) max(|u|, c) is not run here: section
# 18 shows it only sheds contamination as it converges to the quantile, which
# is already in this comparison as alpha = 0.
#
# Output: out/fits-invhuber.rds, out/invhuber-gpd.txt, out/fig-invhuber-gpd.png
# ---------------------------------------------------------------------------
.libPaths(c("/home/user/Rlib-graft", .libPaths()))
source("R/setup.R"); source("R/config.R")
source("R/gpd.R"); source("R/gpd_estimators.R"); source("R/invhuber.R")
source("R/smoothgraft.R"); source("R/graft_fast.R")
suppressMessages({library(distionary); library(distplyr)})
library(parallel)

n <- N_OBS; NC <- detectCores()
Tp <- RETURN_PERIODS; ex <- 1 / Tp; truth_rl <- q_true(1 - ex)
set.seed(4242 + n)                                  # same datasets as scripts/16
datasets <- lapply(seq_len(N_REP), function(i) r_true(n))

W  <- function(p) p^6
WD <- function(p) 6 * p^5
grid0 <- make_level_grid(0, N_PANEL, N_GL)
KS <- c(0.25, 0.5, 1, 2, 4)

big <- r_true(5000)
cat("asymptotic fits (n = 5000) ...\n")
for (k in KS) {
  f <- fit_gpd_onesided(big, grid0, W, k)
  cat(sprintf("  k = %-5.2f  (%.3f, %.3f, %.3f)\n", k, f[1], f[2], f[3]))
}

RL <- list()
RL[["POT-MLE u=0.90"]] <- do.call(rbind, mclapply(datasets,
  function(y) pot_return_levels(y, fit_pot_mle(y, 0.90), ex), mc.cores = NC))
RL[["POT-Lmom u=0.90"]] <- do.call(rbind, mclapply(datasets,
  function(y) pot_return_levels(y, fit_pot_lmom(y, 0.90), ex), mc.cores = NC))

PAR <- list()
for (k in KS) {
  nm <- sprintf("k=%.2f", k)
  cat(sprintf("%s fit ... ", nm)); t0 <- Sys.time()
  PAR[[nm]] <- do.call(rbind, mclapply(datasets,
    function(y) fit_gpd_onesided(y, grid0, W, k), mc.cores = NC))
  cat(sprintf("%.1f min, graft ... ", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  t0 <- Sys.time()
  RL[[paste0(nm, " + graft")]] <- do.call(rbind, mclapply(seq_along(datasets),
    function(i) as.numeric(graft_fast_return_levels(datasets[[i]], PAR[[nm]][i, ],
                                                    W, ex, "gpd", w_deriv = WD)),
    mc.cores = NC))
  cat(sprintf("%.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
saveRDS(list(RL = RL, PAR = PAR, ks = KS, Tp = Tp, truth = truth_rl),
        "out/fits-invhuber.rds")

# pull the pure L1 / L2 / best-elastile rows from scripts/16 for comparison
el <- readRDS("out/fits-gpd-elastile.rds")
for (a in c(0, 0.5, 1)) {
  nm <- sprintf("alpha=%.2f + graft", a)
  RL[[sprintf("[elastile a=%.2f]", a)]] <- el$RL[[nm]]
}

REF <- RL[["POT-MLE u=0.90"]]
summ <- function(r) {
  ok <- complete.cases(r) & complete.cases(REF)
  e  <- sweep(r[ok, , drop = FALSE], 2, truth_rl, "-")
  er <- sweep(REF[ok, , drop = FALSE], 2, truth_rl, "-")
  list(ratio = colMeans(e^2) / colMeans(er^2),
       se = apply(e^2 - er^2, 2, sd) / sqrt(nrow(e)) / colMeans(er^2),
       medae = apply(abs(e), 2, median), win = colMeans(abs(e) < abs(er)),
       bias = colMeans(e), nfail = sum(!ok))
}
S <- lapply(RL, summ)
Tshow <- c(20, 50, 100, 200, 500, 1000)
idx <- sapply(Tshow, function(t) which.min(abs(Tp - t)))

sink("out/invhuber-gpd.txt", split = TRUE)
cat(sprintf("=== One-sided inverted Huber, GPD + graft, n = %d, %d replicates ===\n", n, N_REP))
cat("psi(u) = max(u, -c),  c = k * IQR(y).  w(p) = p^6 fitting and handover.\n")
cat("Ratios against POT-MLE(0.90) hard grafted. Bracketed rows are from\n")
cat("scripts/16 on the same datasets: a=0 pure quantile, a=1 pure expectile.\n\n")
cat("--- median fitted shape ---\n")
print(round(sapply(PAR, function(m) median(m[, 3], na.rm = TRUE)), 3))
for (lab in c("ratio", "se", "medae", "win", "bias")) {
  cat(sprintf("\n--- %s ---\n", switch(lab,
    ratio = "MSE relative to POT-MLE(0.90)", se = "MC standard error of that ratio",
    medae = "median absolute error", bias = "bias",
    win = "fraction closer to the truth than POT-MLE(0.90)")))
  tab <- t(sapply(S, function(s) s[[lab]][idx])); colnames(tab) <- paste0("T=", round(Tp[idx]))
  print(round(tab, 3))
}
cat("\nfailed fits:\n"); print(sapply(S, `[[`, "nfail"))
sink()

png("out/fig-invhuber-gpd.png", width = 1650, height = 640, res = 133)
par(mfrow = c(1, 2), mar = c(4.4, 4.6, 3.4, 1.2))
cols <- colorRampPalette(c("#2980b9", "#7d3c98", "#1e8449"))(length(KS))
Sk <- S[paste0(sprintf("k=%.2f", KS), " + graft")]
plot(Tp, rep(1, length(Tp)), type = "n", log = "xy", ylim = c(0.1, 5),
     xlab = "return period T", ylab = "MSE relative to POT-MLE(0.90)",
     main = sprintf("One-sided inverted Huber + graft (n = %d)", n))
abline(h = 1, col = "grey50", lwd = 1.5)
lines(Tp, S[["[elastile a=1.00]"]]$ratio, col = "#c0392b", lwd = 3, lty = 2)
lines(Tp, S[["[elastile a=0.50]"]]$ratio, col = "#e67e22", lwd = 3, lty = 3)
for (i in seq_along(KS)) lines(Tp, Sk[[i]]$ratio, col = cols[i], lwd = 2.4)
legend("bottomleft", c(sprintf("k = %.2f", KS), "pure L2", "elastile 0.5"),
       col = c(cols, "#c0392b", "#e67e22"), lwd = 2.4,
       lty = c(rep(1, length(KS)), 2, 3), bty = "n", cex = 0.65)
k1 <- which.min(abs(Tp - 1000))
r <- sapply(Sk, function(s) s$ratio[k1])
plot(KS, r, type = "b", pch = 19, lwd = 2.4, col = "#7d3c98", log = "x",
     ylim = range(c(r, S[["[elastile a=1.00]"]]$ratio[k1], 1)),
     xlab = "knot k  (c = k * IQR;  k -> Inf is the expectile)",
     ylab = "MSE relative to POT-MLE(0.90)", main = sprintf("T = %.0f", Tp[k1]))
abline(h = 1, col = "#c0392b", lwd = 2)
abline(h = S[["[elastile a=1.00]"]]$ratio[k1], col = "#c0392b", lwd = 2, lty = 2)
abline(h = S[["[elastile a=0.50]"]]$ratio[k1], col = "#e67e22", lwd = 2, lty = 3)
legend("topright", c("inverted Huber", "POT-MLE(0.90)", "pure L2", "elastile 0.5"),
       col = c("#7d3c98", "#c0392b", "#c0392b", "#e67e22"), lwd = 2.4,
       lty = c(1, 1, 2, 3), bty = "n", cex = 0.7)
dev.off()
cat("\nwrote out/fig-invhuber-gpd.png and out/invhuber-gpd.txt\n")

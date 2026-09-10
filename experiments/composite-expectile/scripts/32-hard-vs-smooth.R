# ---------------------------------------------------------------------------
# 32. Is a hard graft enough, or does the handover have to be smooth?
#
# The composite fits have a poor body by construction -- the weight told them to
# ignore it. Section 9 repairs that with a smooth graft. But a hard graft is far
# simpler: empirical below a threshold, the fitted tail above, rescaled so the
# survival function is continuous. If that already repairs the body, the smooth
# handover is machinery without a purpose.
#
# Three arms on the same 2000 datasets and the same composite fits: ungrafted,
# hard grafted at a sweep of thresholds, and smooth grafted with w = p^6.
#
# Output: out/hard-vs-smooth.txt, out/fig-hard-vs-smooth.png
# ---------------------------------------------------------------------------
.libPaths(c("/home/user/Rlib-graft", .libPaths()))
source("R/setup.R"); source("R/config.R")
source("R/gpd.R"); source("R/gpd_estimators.R"); source("R/invhuber.R")
source("R/hardgraft.R"); source("R/smoothgraft.R"); source("R/graft_fast.R")
suppressMessages({library(distionary); library(distplyr)})
library(parallel)

n <- N_OBS; NC <- detectCores()
Tp <- RETURN_PERIODS; ex <- 1 / Tp; truth_rl <- q_true(1 - ex)
W <- function(p) p^6; WD <- function(p) 6 * p^5
gr <- make_level_grid(0, N_PANEL, N_GL)
set.seed(4242 + n)
dat <- lapply(seq_len(N_REP), function(i) r_true(n))

FITS <- list("composite L2"  = function(y) fit_gpd_composite(y, gr, W, "expectile"),
             "inv Huber k=4" = function(y) fit_gpd_onesided(y, gr, W, 4))
VS <- c(0.50, 0.70, 0.80, 0.90, 0.95)

RL <- list()
RL[["POT-MLE u=0.90"]] <- do.call(rbind, mclapply(dat,
  function(y) pot_return_levels(y, fit_pot_mle(y, 0.90), ex), mc.cores = NC))
PAR <- list()
for (nm in names(FITS)) {
  cat(sprintf("%s fit ... ", nm)); t0 <- Sys.time()
  PAR[[nm]] <- do.call(rbind, mclapply(dat, FITS[[nm]], mc.cores = NC))
  cat(sprintf("%.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  P <- PAR[[nm]]
  RL[[paste0(nm, " | ungrafted")]] <- do.call(rbind, lapply(seq_len(N_REP),
    function(i) if (anyNA(P[i, ])) rep(NA_real_, length(ex)) else
      qgpd(1 - ex, P[i, 1], P[i, 2], P[i, 3])))
  for (v in VS)
    RL[[sprintf("%s | hard v=%.2f", nm, v)]] <- do.call(rbind, mclapply(seq_along(dat),
      function(i) hard_graft_return_levels(dat[[i]], P[i, ], ex, "gpd", v), mc.cores = NC))
  RL[[paste0(nm, " | SMOOTH")]] <- do.call(rbind, mclapply(seq_along(dat),
    function(i) as.numeric(graft_fast_return_levels(dat[[i]], P[i, ], W, ex, "gpd",
                                                    w_deriv = WD)), mc.cores = NC))
}
saveRDS(list(RL = RL, PAR = PAR, Tp = Tp, truth = truth_rl), "out/fits-hard-vs-smooth.rds")

REF <- RL[["POT-MLE u=0.90"]]
idx <- sapply(c(2, 5, 10, 25, 50, 100, 500, 1000), function(t) which.min(abs(Tp - t)))
summ <- function(r) { ok <- complete.cases(r) & complete.cases(REF)
  e <- sweep(r[ok, , drop = FALSE], 2, truth_rl, "-")
  er <- sweep(REF[ok, , drop = FALSE], 2, truth_rl, "-")
  (colMeans(e^2) / colMeans(er^2))[idx] }
S <- lapply(RL, summ)

sink("out/hard-vs-smooth.txt", split = TRUE)
cat(sprintf("=== Hard graft against smooth graft, n = %d, %d replicates ===\n", n, N_REP))
cat("Contaminated truth, MSE relative to POT-MLE(0.90). The hard graft joins an\n")
cat("empirical body to the fitted tail at the v-quantile, rescaled so the survival\n")
cat("function is continuous; the smooth graft blends with w = p^6.\n\n")
cat(sprintf("%-28s", "")); cat(sprintf("%8s", paste0("T=", round(Tp[idx])))); cat("\n")
for (nm in names(S)) { cat(sprintf("%-28s", nm)); cat(sprintf("%8.3f", S[[nm]])); cat("\n") }
sink()

png("out/fig-hard-vs-smooth.png", width = 1700, height = 700, res = 133)
par(mfrow = c(1, 2), mar = c(4.6, 4.8, 3.4, 1.2))
cols <- colorRampPalette(c("#e8b96a", "#a53a2b"))(length(VS))
for (nm in names(FITS)) {
  plot(Tp, rep(1, length(Tp)), type = "n", log = "xy", ylim = c(0.15, 12),
       xlab = "return period T", ylab = "MSE relative to POT-MLE(0.90)", main = nm)
  abline(h = 1, col = "grey55", lwd = 1.5)
  ok <- complete.cases(REF)
  full <- function(r) { o <- complete.cases(r) & ok
    e <- sweep(r[o, , drop = FALSE], 2, truth_rl, "-")
    er <- sweep(REF[o, , drop = FALSE], 2, truth_rl, "-")
    colMeans(e^2) / colMeans(er^2) }
  lines(Tp, full(RL[[paste0(nm, " | ungrafted")]]), col = "#1f5f8b", lwd = 3, lty = 3)
  for (i in seq_along(VS))
    lines(Tp, full(RL[[sprintf("%s | hard v=%.2f", nm, VS[i])]]), col = cols[i], lwd = 2.2)
  lines(Tp, full(RL[[paste0(nm, " | SMOOTH")]]), col = "#197a45", lwd = 3.4)
  legend("topright", c("ungrafted", sprintf("hard, v = %.2f", VS), "smooth graft"),
         col = c("#1f5f8b", cols, "#197a45"), lwd = c(3, rep(2.2, length(VS)), 3.4),
         lty = c(3, rep(1, length(VS)), 1), bty = "n", cex = 0.68)
}
dev.off()
cat("\nwrote out/fig-hard-vs-smooth.png and out/hard-vs-smooth.txt\n")

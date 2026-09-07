# ---------------------------------------------------------------------------
# 24. The same honest test in the GPD setting, which adds a second question.
#
# When the data are iid GPD the whole sample is correctly specified, so
# peaks-over-threshold at ANY threshold is correct too (the GPD is threshold
# stable) -- a lower threshold simply uses more data. That gives two references:
# POT-MLE(0.90), the practitioner default and the reference used throughout
# section 13, and POT-MLE(0.50), which is correctly specified AND more
# efficient.
#
# The second question is the graft. Sections 9 and 13 graft an empirical body
# onto the composite tail because the body is wrong. If the body is exactly
# right, an empirical body is strictly noisier than the correct parametric one,
# so grafting should now HURT. Both are reported to separate the two effects.
#
# Output: out/fits-correct-gpd.rds, out/correct-gpd.txt, out/fig-correct-gpd.png
# ---------------------------------------------------------------------------
.libPaths(c("/home/user/Rlib-graft", .libPaths()))
source("R/setup.R"); source("R/config.R")
source("R/gpd.R"); source("R/gpd_estimators.R"); source("R/invhuber.R")
source("R/smoothgraft.R"); source("R/graft_fast.R")
suppressMessages({library(distionary); library(distplyr)})
library(parallel)

n <- N_OBS; NC <- detectCores()
MU <- 1; SIG <- 0.5; XI <- 0.2                    # the truth, exactly a GPD
Tp <- RETURN_PERIODS; ex <- 1 / Tp
truth_rl <- qgpd(1 - ex, MU, SIG, XI)

set.seed(20260907)
datasets <- lapply(seq_len(N_REP), function(i) qgpd(runif(n), MU, SIG, XI))

W  <- function(p) p^6
WD <- function(p) 6 * p^5
gr <- make_level_grid(0, N_PANEL, N_GL)

rl_par <- function(par) if (anyNA(par)) rep(NA_real_, length(ex)) else
  qgpd(1 - ex, par[1], par[2], par[3])

RL <- list()
for (v in c(0.90, 0.50)) {
  RL[[sprintf("POT-MLE u=%.2f", v)]] <- do.call(rbind, mclapply(datasets,
    function(y) pot_return_levels(y, fit_pot_mle(y, v), ex), mc.cores = NC))
}
RL[["POT-Lmom u=0.90"]] <- do.call(rbind, mclapply(datasets,
  function(y) pot_return_levels(y, fit_pot_lmom(y, 0.90), ex), mc.cores = NC))

FITS <- list(
  "composite L1"   = function(y) fit_gpd_composite(y, gr, W, "quantile"),
  "composite L2"   = function(y) fit_gpd_composite(y, gr, W, "expectile"),
  "elastile a=0.5" = function(y) fit_gpd_elastile(y, gr, W, 0.5),
  "inv Huber k=4"  = function(y) fit_gpd_onesided(y, gr, W, 4)
)
PAR <- list()
for (nm in names(FITS)) {
  cat(sprintf("%-15s fit ... ", nm)); t0 <- Sys.time()
  PAR[[nm]] <- do.call(rbind, mclapply(datasets, FITS[[nm]], mc.cores = NC))
  RL[[nm]] <- do.call(rbind, lapply(seq_len(N_REP), function(i) rl_par(PAR[[nm]][i, ])))
  cat(sprintf("%.1f min, graft ... ", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  t0 <- Sys.time()
  RL[[paste0(nm, " + graft")]] <- do.call(rbind, mclapply(seq_along(datasets),
    function(i) as.numeric(graft_fast_return_levels(datasets[[i]], PAR[[nm]][i, ],
                                                    W, ex, "gpd", w_deriv = WD)),
    mc.cores = NC))
  cat(sprintf("%.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
saveRDS(list(PAR = PAR, RL = RL, Tp = Tp, truth = truth_rl,
             pars = c(MU, SIG, XI)), "out/fits-correct-gpd.rds")

REF <- RL[["POT-MLE u=0.90"]]
summ <- function(r) {
  ok <- complete.cases(r) & complete.cases(REF)
  e <- sweep(r[ok, , drop = FALSE], 2, truth_rl, "-")
  er <- sweep(REF[ok, , drop = FALSE], 2, truth_rl, "-")
  list(ratio = colMeans(e^2) / colMeans(er^2),
       medae = apply(abs(e), 2, median), win = colMeans(abs(e) < abs(er)),
       bias = colMeans(e), nfail = sum(!ok))
}
S <- lapply(RL, summ)
idx <- sapply(c(2, 10, 50, 100, 200, 500, 1000), function(t) which.min(abs(Tp - t)))

sink("out/correct-gpd.txt", split = TRUE)
cat(sprintf("=== Correctly specified GPD: iid GPD(%g, %g, %g), n = %d, %d reps ===\n",
            MU, SIG, XI, n, N_REP))
cat("Ratios against POT-MLE(0.90), the reference used in section 13. The GPD is\n")
cat("threshold stable, so POT-MLE(0.50) is also correct and uses more data.\n\n")
cat("--- median fitted parameters (truth 1, 0.5, 0.2) ---\n")
print(round(t(sapply(PAR, function(m) apply(m, 2, median, na.rm = TRUE))), 3))
for (lab in c("ratio", "medae", "win", "bias")) {
  cat(sprintf("\n--- %s ---\n", switch(lab,
    ratio = "MSE relative to POT-MLE(0.90)", medae = "median absolute error",
    win = "fraction closer to the truth than POT-MLE(0.90)", bias = "bias")))
  tab <- t(sapply(S, function(s) s[[lab]][idx])); colnames(tab) <- paste0("T=", round(Tp[idx]))
  print(round(tab, 3))
}
cat("\nfailed fits:\n"); print(sapply(S, `[[`, "nfail"))
sink()

png("out/fig-correct-gpd.png", width = 1650, height = 620, res = 133)
par(mfrow = c(1, 2), mar = c(4.6, 4.6, 3.4, 1.2))
nmk <- names(FITS); cols <- c("#1f5f8b", "#197a45", "#7d3c98", "#a53a2b")
for (g in c(FALSE, TRUE)) {
  plot(Tp, rep(1, length(Tp)), type = "n", log = "x", ylim = c(0.5, 3),
       xlab = "return period T", ylab = "MSE relative to POT-MLE(0.90)",
       main = if (g) "grafted onto an empirical body" else "ungrafted")
  abline(h = 1, col = "grey40", lwd = 2)
  lines(Tp, S[["POT-MLE u=0.50"]]$ratio, col = "black", lwd = 3, lty = 2)
  lines(Tp, S[["POT-Lmom u=0.90"]]$ratio, col = "#e67e22", lwd = 2.4, lty = 3)
  for (i in seq_along(nmk))
    lines(Tp, S[[if (g) paste0(nmk[i], " + graft") else nmk[i]]]$ratio,
          col = cols[i], lwd = 2.4)
  legend("topleft", c(nmk, "POT-MLE(0.50)", "POT-Lmom(0.90)"),
         col = c(cols, "black", "#e67e22"), lwd = 2.4,
         lty = c(rep(1, 4), 2, 3), bty = "n", cex = 0.65)
}
dev.off()
cat("\nwrote out/fig-correct-gpd.png and out/correct-gpd.txt\n")

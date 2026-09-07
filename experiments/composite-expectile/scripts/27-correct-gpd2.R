# ---------------------------------------------------------------------------
# 27. The correctly-specified GPD test, done properly: truly just a GPD.
#
# scripts/24-26 fitted three parameters (mu, sigma, xi) to data drawn from
# GPD(1, 0.5, 0.2). Two things were wrong with that.
#
#   1. A three-parameter GPD with unknown threshold is non-regular: the
#      likelihood is monotone in mu up to min(y), so mu_hat is a boundary
#      estimate. Worse for the comparison, the "oracle" reference was given the
#      true threshold while the composite estimators had to estimate it, so they
#      were paying for a parameter the reference got free.
#
#   2. A GPD used as a tail model does not have a location to estimate. It
#      starts at the threshold, and only the scale and shape are unknown.
#
# So: draw iid GPD(0, sigma, xi) and fit only (sigma, xi), for every estimator
# including the references. Everything now has the same two unknowns and the
# problem is regular, so the MLE is genuinely efficient.
#
# Also included is the empirical distribution on its own, to settle what the
# short-return-period numbers are actually measuring.
#
# ON RETURN PERIODS. Here the GPD IS the whole distribution, so an exceedance
# probability of 1/T is a T-year return level and the labels are self
# consistent. In a peaks-over-threshold deployment the GPD sits above a
# threshold with exceedance rate zeta, and the T-year level solves
# zeta * S_gpd(x) = 1/T, i.e. S_gpd(x) = 1/(zeta T). So a column labelled T here
# is the (T / zeta)-year level in that setting: with zeta = 0.1, the column
# labelled T = 100 is the 1000-year return level. The MSE ratios are unchanged
# -- only the label moves -- but it moves by a factor of ten, so the
# practically relevant columns are the middle of this table, not the right edge.
#
# Output: out/fits-correct-gpd2.rds, out/correct-gpd2.txt, out/fig-correct-gpd2.png
# ---------------------------------------------------------------------------
.libPaths(c("/home/user/Rlib-graft", .libPaths()))
source("R/setup.R"); source("R/config.R")
source("R/gpd.R"); source("R/gpd_estimators.R"); source("R/invhuber.R"); source("R/gpd2.R")
source("R/smoothgraft.R"); source("R/graft_fast.R")
suppressMessages({library(distionary); library(distplyr)})
library(parallel)

n <- N_OBS; NC <- detectCores()
SIG <- 0.5; XI <- 0.2                              # truth: GPD(0, 0.5, 0.2)
Tp <- RETURN_PERIODS; ex <- 1 / Tp
truth_rl <- qgpd(1 - ex, 0, SIG, XI)
ZETA <- 0.1                                        # for the POT relabelling only

set.seed(20260908)
datasets <- lapply(seq_len(N_REP), function(i) qgpd(runif(n), 0, SIG, XI))

W <- function(p) p^6; WD <- function(p) 6 * p^5
gr <- make_level_grid(0, N_PANEL, N_GL)
rl2 <- function(par) if (anyNA(par)) rep(NA_real_, length(ex)) else
  qgpd(1 - ex, 0, par[1], par[2])

RL <- list()
RL[["GPD MLE"]]       <- do.call(rbind, mclapply(datasets, function(y) rl2(gpd2_mle(y)), mc.cores = NC))
RL[["GPD L-moments"]] <- do.call(rbind, mclapply(datasets, function(y) rl2(gpd2_lmom(y)), mc.cores = NC))
# the empirical distribution alone: it cannot extrapolate past the sample max
RL[["empirical only"]] <- do.call(rbind, mclapply(datasets,
  function(y) as.numeric(stats::quantile(y, 1 - ex, type = 7)), mc.cores = NC))

FITS <- list(
  "composite L1"   = function(y) fit_gpd2_composite(y, gr, W, "quantile"),
  "composite L2"   = function(y) fit_gpd2_composite(y, gr, W, "expectile"),
  "elastile a=0.5" = function(y) fit_gpd2_elastile(y, gr, W, 0.5),
  "inv Huber k=4"  = function(y) fit_gpd2_onesided(y, gr, W, 4)
)
PAR <- list()
for (nm in names(FITS)) {
  cat(sprintf("%-15s fit ... ", nm)); t0 <- Sys.time()
  PAR[[nm]] <- do.call(rbind, mclapply(datasets, FITS[[nm]], mc.cores = NC))
  RL[[nm]] <- do.call(rbind, lapply(seq_len(N_REP), function(i) rl2(PAR[[nm]][i, ])))
  cat(sprintf("%.1f min, graft ... ", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  t0 <- Sys.time()
  RL[[paste0(nm, " + graft")]] <- do.call(rbind, mclapply(seq_along(datasets),
    function(i) as.numeric(graft_fast_return_levels(datasets[[i]],
      c(0, PAR[[nm]][i, 1], PAR[[nm]][i, 2]), W, ex, "gpd", w_deriv = WD)),
    mc.cores = NC))
  cat(sprintf("%.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
saveRDS(list(RL = RL, PAR = PAR, Tp = Tp, truth = truth_rl, pars = c(SIG, XI)),
        "out/fits-correct-gpd2.rds")

REF <- RL[["GPD MLE"]]
summ <- function(r) {
  ok <- complete.cases(r) & complete.cases(REF)
  e  <- sweep(r[ok, , drop = FALSE], 2, truth_rl, "-")
  er <- sweep(REF[ok, , drop = FALSE], 2, truth_rl, "-")
  list(ratio = colMeans(e^2) / colMeans(er^2), bias = colMeans(e),
       win = colMeans(abs(e) < abs(er)), nfail = sum(!ok))
}
S <- lapply(RL, summ)
idx <- sapply(c(2, 5, 10, 25, 50, 100, 500, 1000), function(t) which.min(abs(Tp - t)))

sink("out/correct-gpd2.txt", split = TRUE)
cat(sprintf("=== Truly just a GPD: iid GPD(0, %g, %g), mu KNOWN, n = %d, %d reps ===\n",
            SIG, XI, n, N_REP))
cat("Only scale and shape are estimated, by every method including the MLE, so\n")
cat("the problem is regular and the MLE is efficient. Ratios are against it.\n\n")
cat("Column labels are the GPD's OWN return period. Attached above a threshold\n")
cat(sprintf("with exceedance rate zeta = %g, the column labelled T is the %g-year\n", ZETA, 1/ZETA))
cat("level times T -- so T = 100 here is the 1000-year level in that setting.\n\n")
cat("--- median fitted (sigma, xi); truth (0.5, 0.2) ---\n")
print(round(t(sapply(PAR, function(m) apply(m, 2, median, na.rm = TRUE))), 3))
for (lab in c("ratio", "win", "bias")) {
  cat(sprintf("\n--- %s ---\n", switch(lab, ratio = "MSE relative to the GPD MLE",
    win = "fraction closer to the truth than the MLE", bias = "bias")))
  tab <- t(sapply(S, function(s) s[[lab]][idx]))
  colnames(tab) <- paste0("T=", round(Tp[idx])); print(round(tab, 3))
}
cat(sprintf("\n(equivalently, with zeta = %g those columns are T = %s)\n",
            ZETA, paste(round(Tp[idx] / ZETA), collapse = ", ")))
cat("\nfailed fits:\n"); print(sapply(S, `[[`, "nfail"))
sink()

png("out/fig-correct-gpd2.png", width = 1650, height = 620, res = 133)
par(mfrow = c(1, 2), mar = c(4.6, 4.6, 3.4, 1.2))
cols <- c("#1f5f8b", "#197a45", "#7d3c98", "#a53a2b")
nmk <- names(FITS)
for (g in c(FALSE, TRUE)) {
  plot(Tp, rep(1, length(Tp)), type = "n", log = "x", ylim = c(0.6, 4),
       xlab = "GPD's own return period T   (zeta = 0.1: multiply by 10)",
       ylab = "MSE relative to the GPD MLE",
       main = if (g) "grafted onto an empirical body" else "ungrafted")
  abline(h = 1, col = "grey40", lwd = 2)
  lines(Tp, S[["GPD L-moments"]]$ratio, col = "#e67e22", lwd = 3, lty = 2)
  for (i in seq_along(nmk))
    lines(Tp, S[[if (g) paste0(nmk[i], " + graft") else nmk[i]]]$ratio,
          col = cols[i], lwd = 2.4)
  if (g) lines(Tp, S[["empirical only"]]$ratio, col = "grey30", lwd = 2.4, lty = 3)
  legend("topleft", c(nmk, "L-moments", if (g) "empirical only"),
         col = c(cols, "#e67e22", if (g) "grey30"), lwd = 2.4,
         lty = c(rep(1, 4), 2, if (g) 3), bty = "n", cex = 0.62)
}
dev.off()
cat("\nwrote out/fig-correct-gpd2.png and out/correct-gpd2.txt\n")

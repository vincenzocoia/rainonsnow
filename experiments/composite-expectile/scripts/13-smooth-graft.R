# ---------------------------------------------------------------------------
# Does grafting the composite fit onto an empirical body fix the ruined body?
#
# The composite estimators buy their tail fit by abandoning the body: at
# p0 = 0.95 the 2-year return level has a mean squared error ten thousand times
# the MLE's. The proposal is to keep the composite fit only as the TAIL of a
# smooth graft (Coia's hazard-mixture construction, probaverse/distplyr branch
# claude/distplyr-smooth-graft-spond9) whose body is the empirical distribution,
# handing over with the same weight function used in the fitting criterion.
#
# Three things are separated here:
#   - the ungrafted composite fits, as a reference;
#   - grafts that hand over with the SAME weight used for fitting;
#   - grafts that hand over on a fixed [0.90, 0.98] window regardless of the
#     fitting weight, since the two weights do different jobs;
# and, as the control that decides whether any of this is about the composite
# estimator at all, grafts of the MLE and L-moments fits onto the same body.
#
# Output: out/fits-graftfit.rds, out/smooth-graft.txt, out/fig-smooth-graft.png
# ---------------------------------------------------------------------------
.libPaths(c("/home/user/Rlib-graft", .libPaths()))   # distionary(dev) + distplyr
source("R/setup.R")
source("R/config.R")
source("R/smoothgraft.R")
suppressMessages({library(distionary); library(distplyr)})
library(parallel)

n <- N_OBS
Tp <- RETURN_PERIODS; ex <- 1 / Tp; truth_rl <- q_true(1 - ex)
set.seed(4242 + n)                       # the datasets used in scripts/01 and 08
datasets <- lapply(seq_len(N_REP), function(i) r_true(n))

base <- readRDS(sprintf("out/fits-n%d.rds", n))$fits
tf   <- readRDS("out/fits-tailfocus.rds")$fits
PAR <- list(mle = base$mle, lmom = base$lmom,
            `cqe|0.95` = base$cqe, `cee|0.95` = base$cee,
            `cqe|0.50` = tf[["quantile|0.50"]], `cee|0.50` = tf[["expectile|0.50"]])

# Grafted variants: which fit supplies the tail, and where the handover runs.
GRAFTS <- list(
  `graft cee|0.95, handover 0.95-0.98` = list(par = "cee|0.95", p0 = 0.95, p1 = 0.98),
  `graft cqe|0.95, handover 0.95-0.98` = list(par = "cqe|0.95", p0 = 0.95, p1 = 0.98),
  `graft cee|0.50, handover 0.50-0.53` = list(par = "cee|0.50", p0 = 0.50, p1 = 0.53),
  `graft cee|0.50, handover 0.90-0.98` = list(par = "cee|0.50", p0 = 0.90, p1 = 0.98),
  `graft cqe|0.50, handover 0.90-0.98` = list(par = "cqe|0.50", p0 = 0.90, p1 = 0.98),
  `graft MLE,      handover 0.90-0.98` = list(par = "mle",      p0 = 0.90, p1 = 0.98),
  `graft Lmom,     handover 0.90-0.98` = list(par = "lmom",     p0 = 0.90, p1 = 0.98)
)

RL <- list()
# Parametric fits, ungrafted.
for (nm in names(PAR))
  RL[[nm]] <- t(apply(PAR[[nm]], 1, function(th)
    if (any(is.na(th))) rep(NA_real_, length(Tp)) else qgev(1 - ex, th[1], th[2], th[3])))
# The empirical distribution on its own, as the body baseline.
cat("empirical ... ")
RL[["empirical"]] <- do.call(rbind, mclapply(datasets, function(y)
  as.numeric(eval_quantile(dst_empirical(y), at = 1 - ex)), mc.cores = detectCores()))
cat("done\n")
# The grafts.
for (nm in names(GRAFTS)) {
  g <- GRAFTS[[nm]]; P <- PAR[[g$par]]
  cat(sprintf("%-38s ... ", nm)); t0 <- Sys.time()
  RL[[nm]] <- do.call(rbind, mclapply(seq_along(datasets), function(i)
    as.numeric(graft_return_levels(datasets[[i]], P[i, ], g$p0, g$p1, ex)),
    mc.cores = detectCores()))
  cat(sprintf("%.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
saveRDS(list(RL = RL, Tp = Tp, truth = truth_rl, n = n), "out/fits-graftfit.rds")

M <- RL$mle
summ <- function(r) {
  ok <- complete.cases(r) & complete.cases(M)
  r <- r[ok, , drop = FALSE]; m <- M[ok, , drop = FALSE]
  e <- sweep(r, 2, truth_rl, "-"); em <- sweep(m, 2, truth_rl, "-")
  list(bias = colMeans(e), sd = apply(r, 2, sd), mse = colMeans(e^2),
       medae = apply(abs(e), 2, median),
       mse_ratio = colMeans(e^2) / colMeans(em^2),
       mse_se = apply(e^2 - em^2, 2, sd) / sqrt(nrow(e)) / colMeans(em^2),
       win = colMeans(abs(e) < abs(em)), nfail = sum(!ok))
}
S <- lapply(RL, summ)
Tshow <- c(2, 5, 10, 20, 50, 100, 200, 500, 1000)
idx <- sapply(Tshow, function(t) which.min(abs(Tp - t)))

sink("out/smooth-graft.txt", split = TRUE)
cat(sprintf("=== Smooth graft: empirical body + composite tail, n = %d, %d replicates ===\n\n",
            n, N_REP))
for (lab in c("mse_ratio", "medae", "bias", "sd", "win")) {
  cat(sprintf("\n--- %s ---\n", switch(lab,
    mse_ratio = "MSE relative to GEV MLE", medae = "median absolute error",
    bias = "bias", sd = "sd", win = "fraction of datasets closer to the truth than the MLE")))
  tab <- t(sapply(S, function(s) s[[lab]][idx]))
  colnames(tab) <- paste0("T=", round(Tp[idx])); print(round(tab, 3))
}
cat("\n--- Monte Carlo standard error of each MSE ratio ---\n")
tab <- t(sapply(S, function(s) s$mse_se[idx])); colnames(tab) <- paste0("T=", round(Tp[idx]))
print(round(tab, 3))
cat("\nfailed fits per variant:\n"); print(sapply(S, `[[`, "nfail"))
sink()

png("out/fig-smooth-graft.png", width = 1650, height = 660, res = 133)
par(mfrow = c(1, 2), mar = c(4.4, 4.6, 3.4, 1.2))
show <- c("cee|0.95", "graft cee|0.95, handover 0.95-0.98",
          "cee|0.50", "graft cee|0.50, handover 0.90-0.98",
          "empirical", "graft MLE,      handover 0.90-0.98", "lmom")
cols <- c("#1e8449", "#1e8449", "#2980b9", "#2980b9", "#7f8c8d", "#c0392b", "#e67e22")
ltys <- c(2, 1, 2, 1, 3, 1, 4)
plot(Tp, rep(1, length(Tp)), type = "n", log = "xy", ylim = c(0.2, 3e4),
     xlab = "return period T (years)", ylab = "MSE relative to GEV MLE",
     main = sprintf("Grafting the composite fit onto\nan empirical body (n = %d)", n))
abline(h = 1, col = "grey50", lwd = 1.5)
for (i in seq_along(show)) lines(Tp, S[[show[i]]]$mse_ratio, col = cols[i], lwd = 2.6, lty = ltys[i])
legend("topright", show, col = cols, lwd = 2.6, lty = ltys, bty = "n", cex = 0.62)
plot(Tp, rep(1, length(Tp)), type = "n", log = "x", ylim = c(0.3, 1.6),
     xlab = "return period T (years)", ylab = "median |error| relative to GEV MLE",
     main = "Median absolute error\n(robust to the heavy MSE tail)")
abline(h = 1, col = "grey50", lwd = 1.5)
for (i in seq_along(show)) lines(Tp, S[[show[i]]]$medae / S$mle$medae, col = cols[i],
                                 lwd = 2.6, lty = ltys[i])
dev.off()
cat("\nwrote out/fig-smooth-graft.png and out/smooth-graft.txt\n")

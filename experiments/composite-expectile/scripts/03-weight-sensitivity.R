# ---------------------------------------------------------------------------
# How much does the answer depend on where the weight switches on?
#
# p0 is the one real design knob. Moving it up buys less contamination bias at
# the cost of less data, so it dials bias against variance. This re-runs both
# composite estimators at three settings, at the primary sample size.
#
# Output: out/fits-p0.rds, out/fig-p0.png, out/p0-sensitivity.txt
# ---------------------------------------------------------------------------
source("R/setup.R")
source("R/config.R")
library(parallel)

P0_SET <- c(0.90, 0.95, 0.98)
set.seed(4242 + N_OBS)                       # same datasets as scripts/01
datasets <- lapply(seq_len(N_REP), function(i) r_true(N_OBS))

out <- list()
for (p0 in P0_SET) {
  g <- make_level_grid(p0, N_PANEL, N_GL); wf <- make_weight(p0, min(p0 + 0.03, 0.999))
  for (ty in c("quantile", "expectile")) {
    cat(sprintf("p0 = %.2f, %s ... ", p0, ty)); t0 <- Sys.time()
    m <- do.call(rbind, mclapply(datasets, function(y)
      fit_composite(y, g, wf, ty, start = fit_lmom(y)), mc.cores = detectCores()))
    cat(sprintf("%.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
    out[[sprintf("%s_p0_%.2f", ty, p0)]] <- m
  }
}
saveRDS(list(p0 = P0_SET, fits = out, n = N_OBS), "out/fits-p0.rds")

# --- summarise --------------------------------------------------------------
Tp <- RETURN_PERIODS; truth_rl <- q_true(1 - 1 / Tp)
summ <- function(par_mat) {
  rl <- t(apply(par_mat, 1, function(th) if (any(is.na(th))) rep(NA_real_, length(Tp))
                else qgev(1 - 1 / Tp, th[1], th[2], th[3])))
  rl <- rl[complete.cases(rl), , drop = FALSE]
  err <- sweep(rl, 2, truth_rl, "-")
  list(bias = colMeans(err), var = apply(rl, 2, var), mse = colMeans(err^2))
}
S <- lapply(out, summ)
base <- readRDS(sprintf("out/fits-n%d.rds", N_OBS))
mle_mse <- summ(base$fits$mle)$mse

sink("out/p0-sensitivity.txt", split = TRUE)
cat("=== MSE relative to GEV MLE, by where the weight switches on ===\n")
cat(sprintf("n = %d, %d replicates. Values below 1 beat the MLE.\n\n", N_OBS, N_REP))
Tshow <- c(10, 20, 50, 100, 200, 500, 1000)
idx <- sapply(Tshow, function(t) which.min(abs(Tp - t)))
tab <- t(sapply(names(S), function(k) (S[[k]]$mse / mle_mse)[idx]))
colnames(tab) <- paste0("T=", round(Tp[idx]))
print(round(tab, 3))
cat("\n=== bias in the return level ===\n")
tb <- t(sapply(names(S), function(k) S[[k]]$bias[idx])); colnames(tb) <- colnames(tab)
print(round(tb, 3))
cat("\n=== sd of the return level ===\n")
tv <- t(sapply(names(S), function(k) sqrt(S[[k]]$var[idx]))); colnames(tv) <- colnames(tab)
print(round(tv, 3))
sink()

png("out/fig-p0.png", width = 1600, height = 640, res = 135)
par(mfrow = c(1, 2), mar = c(4.4, 4.6, 3.2, 1.2))
cols <- c("#a9cce3", "#2980b9", "#154360")
for (ty in c("quantile", "expectile")) {
  plot(Tp, rep(1, length(Tp)), type = "n", log = "x", ylim = c(0, 3),
       xlab = "return period T", ylab = "MSE relative to GEV MLE",
       main = sprintf("Composite %s (%s)", ty, if (ty == "quantile") "L1" else "L2"))
  abline(h = 1, col = "grey60")
  for (i in seq_along(P0_SET))
    lines(Tp, S[[sprintf("%s_p0_%.2f", ty, P0_SET[i])]]$mse / mle_mse,
          col = cols[i], lwd = 2.6)
  legend("topright", sprintf("p0 = %.2f", P0_SET), col = cols, lwd = 2.6, bty = "n")
}
dev.off()
cat("\nwrote out/fig-p0.png and out/p0-sensitivity.txt\n")

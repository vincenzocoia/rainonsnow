# ---------------------------------------------------------------------------
# Where the weight function should live.
#
# Two dials. The lower edge p0 decides how much contaminated body the estimator
# is allowed to see. The upper edge decides how far past the data it
# extrapolates: with the weight left at one all the way to p = 1, a large share
# of its mass sits above the largest observation, where the empirical functional
# has saturated and the loss can only pull the fitted tail down (see
# scripts/04-shape-diagnostic.R). Capping the weight at the level of the largest
# observation, 1 - 1/n, removes that.
#
# Output: out/fits-weight.rds, out/fig-weight-design.png, out/weight-design.txt
# ---------------------------------------------------------------------------
source("R/setup.R")
source("R/config.R")
library(parallel)

n <- N_OBS
DESIGNS <- list(
  `p0=0.90, uncapped` = list(p0 = 0.90, p1 = 0.93, p2 = 1,          p3 = 1),
  `p0=0.95, uncapped` = list(p0 = 0.95, p1 = 0.98, p2 = 1,          p3 = 1),
  `p0=0.98, uncapped` = list(p0 = 0.98, p1 = 0.99, p2 = 1,          p3 = 1),
  `p0=0.90, capped`   = list(p0 = 0.90, p1 = 0.93, p2 = 1 - 2 / n,  p3 = 1 - 1 / n),
  `p0=0.95, capped`   = list(p0 = 0.95, p1 = 0.97, p2 = 1 - 2 / n,  p3 = 1 - 1 / n)
)

set.seed(4242 + n)                            # the same datasets as scripts/01
datasets <- lapply(seq_len(N_REP), function(i) r_true(n))

fits <- list()
for (nm in names(DESIGNS)) {
  d <- DESIGNS[[nm]]
  wf <- make_weight(d$p0, d$p1, d$p2, d$p3)
  g  <- make_level_grid(d$p0, N_PANEL, N_GL)
  for (ty in c("quantile", "expectile")) {
    cat(sprintf("%-20s %-9s ... ", nm, ty)); t0 <- Sys.time()
    fits[[paste(nm, ty, sep = " | ")]] <- do.call(rbind, mclapply(datasets,
      function(y) fit_composite(y, g, wf, ty, start = fit_lmom(y)),
      mc.cores = detectCores()))
    cat(sprintf("%.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  }
}
saveRDS(list(designs = DESIGNS, fits = fits, n = n), "out/fits-weight.rds")

Tp <- RETURN_PERIODS; truth_rl <- q_true(1 - 1 / Tp)
summ <- function(m) {
  rl <- t(apply(m, 1, function(th) if (any(is.na(th))) rep(NA_real_, length(Tp))
                else qgev(1 - 1 / Tp, th[1], th[2], th[3])))
  rl <- rl[complete.cases(rl), , drop = FALSE]
  err <- sweep(rl, 2, truth_rl, "-")
  list(bias = colMeans(err), sd = apply(rl, 2, sd), mse = colMeans(err^2),
       xi_med = median(m[, 3], na.rm = TRUE),
       at_bound = mean(m[, 3] <= XI_LO + 1e-6, na.rm = TRUE))
}
S <- lapply(fits, summ)
base <- readRDS(sprintf("out/fits-n%d.rds", n))
mle <- summ(base$fits$mle)

sink("out/weight-design.txt", split = TRUE)
cat(sprintf("=== Weight design, n = %d, %d replicates ===\n", n, N_REP))
cat("'capped' means the weight falls back to zero by p = 1 - 1/n, the level of\n")
cat("the largest observation, instead of staying at one up to p = 1.\n\n")
Tshow <- c(10, 20, 50, 100, 200, 500, 1000)
idx <- sapply(Tshow, function(t) which.min(abs(Tp - t)))
cat("--- median fitted shape (true xi = 0.2) and fraction stuck on the shape bound ---\n")
print(round(data.frame(xi_median = sapply(S, `[[`, "xi_med"),
                       frac_at_bound = sapply(S, `[[`, "at_bound")), 3))
cat("\n--- MSE of the return level, relative to GEV MLE (below 1 beats the MLE) ---\n")
tab <- t(sapply(S, function(s) (s$mse / mle$mse)[idx]))
colnames(tab) <- paste0("T=", round(Tp[idx])); print(round(tab, 3))
cat("\n--- bias ---\n")
tb <- t(sapply(S, function(s) s$bias[idx])); colnames(tb) <- colnames(tab); print(round(tb, 3))
cat("\n--- sd ---\n")
tv <- t(sapply(S, function(s) s$sd[idx])); colnames(tv) <- colnames(tab); print(round(tv, 3))
cat("\n--- GEV MLE, for reference ---\n")
print(round(data.frame(bias = mle$bias[idx], sd = mle$sd[idx], mse = mle$mse[idx],
                       row.names = colnames(tab)), 3))
sink()

png("out/fig-weight-design.png", width = 1600, height = 660, res = 135)
par(mfrow = c(1, 2), mar = c(4.4, 4.6, 3.4, 1.2))
cols <- c("#a9cce3", "#2980b9", "#154360", "#f5b041", "#b9770e")
for (ty in c("quantile", "expectile")) {
  plot(Tp, rep(1, length(Tp)), type = "n", log = "xy", ylim = c(0.2, 60),
       xlab = "return period T (years)", ylab = "MSE relative to GEV MLE",
       main = sprintf("Composite %s (%s), n = %d", ty,
                      if (ty == "quantile") "L1" else "L2", n))
  abline(h = 1, col = "grey50", lwd = 1.5)
  for (i in seq_along(DESIGNS))
    lines(Tp, S[[paste(names(DESIGNS)[i], ty, sep = " | ")]]$mse / mle$mse,
          col = cols[i], lwd = 2.6, lty = if (i <= 3) 1 else 2)
  legend("topright", names(DESIGNS), col = cols, lwd = 2.6,
         lty = c(1, 1, 1, 2, 2), bty = "n", cex = 0.72)
}
dev.off()
cat("\nwrote out/fig-weight-design.png and out/weight-design.txt\n")

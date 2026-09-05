# ---------------------------------------------------------------------------
# The headline figures: the four estimators as specified, plus the two
# composite estimators with the weight capped so it does not extend past the
# data. Requires scripts/01 and scripts/03 to have run.
#
# Output: out/fig-headline-return-level.png, out/fig-headline-mse.png
# ---------------------------------------------------------------------------
source("R/setup.R")
source("R/config.R")

Tp <- RETURN_PERIODS; pT <- 1 - 1 / Tp; truth_rl <- q_true(pT)
rls <- function(m) {
  r <- t(apply(m, 1, function(th) if (any(is.na(th))) rep(NA_real_, length(Tp))
               else qgev(pT, th[1], th[2], th[3])))
  r[complete.cases(r), , drop = FALSE]
}
summ <- function(m) { r <- rls(m); e <- sweep(r, 2, truth_rl, "-")
  list(median = apply(r, 2, median), mse = colMeans(e^2), bias = colMeans(e),
       sd = apply(r, 2, sd),
       lo = apply(r, 2, quantile, .025), hi = apply(r, 2, quantile, .975)) }

base <- readRDS(sprintf("out/fits-n%d.rds", N_OBS))$fits
wt   <- readRDS("out/fits-weight.rds")$fits
P <- list(
  list(k = "GEV MLE",                        m = base$mle,  col = "#c0392b"),
  list(k = "GEV L-moments",                  m = base$lmom, col = "#e67e22"),
  list(k = "Composite quantile (L1)",        m = base$cqe,  col = "#2980b9"),
  list(k = "Composite expectile (L2)",       m = base$cee,  col = "#1e8449"),
  list(k = "Composite quantile (L1), capped",
       m = wt[["p0=0.95, capped | quantile"]],  col = "#5dade2"),
  list(k = "Composite expectile (L2), capped",
       m = wt[["p0=0.95, capped | expectile"]], col = "#52be80"))
S <- lapply(P, function(x) summ(x$m))

png("out/fig-headline-return-level.png", width = 1750, height = 1080, res = 132)
par(mfrow = c(2, 3), mar = c(4.2, 4.3, 3.0, 1.0))
for (i in seq_along(P)) {
  plot(Tp, truth_rl, type = "n", log = "x", ylim = c(0, 32),
       xlab = "return period T (years)", ylab = "return level",
       main = sprintf("%s\n(n = %d)", P[[i]]$k, N_OBS), cex.main = 0.95)
  polygon(c(Tp, rev(Tp)), c(S[[i]]$lo, rev(S[[i]]$hi)),
          col = adjustcolor(P[[i]]$col, 0.20), border = NA)
  lines(Tp, S[[i]]$lo, col = P[[i]]$col, lty = 3); lines(Tp, S[[i]]$hi, col = P[[i]]$col, lty = 3)
  lines(Tp, S[[i]]$median, col = P[[i]]$col, lwd = 2.6)
  lines(Tp, truth_rl, col = "black", lwd = 2.6)
  if (i == 1) legend("topleft", c("truth", "median estimate", "central 95%"),
    col = c("black", P[[1]]$col, adjustcolor(P[[1]]$col, .35)), lwd = c(2.6, 2.6, 8),
    bty = "n", cex = 0.8)
}
dev.off()

png("out/fig-headline-mse.png", width = 1600, height = 700, res = 133)
par(mfrow = c(1, 2), mar = c(4.4, 4.6, 3.2, 1.0))
mse <- sapply(S, `[[`, "mse")
plot(Tp, mse[, 1], type = "n", log = "xy", ylim = range(mse),
     xlab = "return period T (years)", ylab = "MSE of estimated return level",
     main = sprintf("MSE over return period (n = %d)", N_OBS))
for (i in seq_along(P)) lines(Tp, S[[i]]$mse, col = P[[i]]$col, lwd = 2.6,
                              lty = if (i > 4) 5 else 1)
legend("bottomright", sapply(P, `[[`, "k"), col = sapply(P, `[[`, "col"),
       lwd = 2.6, lty = c(1,1,1,1,5,5), bty = "n", cex = 0.72)
plot(Tp, mse[, 1] / mse[, 1], type = "n", log = "x", ylim = c(0, 3),
     xlab = "return period T (years)", ylab = "MSE relative to GEV MLE",
     main = "MSE ratio (below 1 beats the MLE)")
abline(h = 1, col = "grey55", lwd = 1.5)
for (i in seq_along(P)) lines(Tp, S[[i]]$mse / mse[, 1], col = P[[i]]$col, lwd = 2.6,
                              lty = if (i > 4) 5 else 1)
legend("topright", sapply(P, `[[`, "k"), col = sapply(P, `[[`, "col"),
       lwd = 2.6, lty = c(1,1,1,1,5,5), bty = "n", cex = 0.72)
dev.off()

sink("out/headline.txt", split = TRUE)
cat(sprintf("=== Headline comparison, n = %d, %d replicates ===\n\n", N_OBS, N_REP))
Tshow <- c(2, 10, 20, 50, 100, 200, 500, 1000)
idx <- sapply(Tshow, function(t) which.min(abs(Tp - t)))
nm <- sapply(P, `[[`, "k")
for (st in c("bias", "sd", "mse")) {
  cat("\n", st, "\n", sep = "")
  tab <- t(sapply(S, function(s) s[[st]][idx])); rownames(tab) <- nm
  colnames(tab) <- paste0("T=", round(Tp[idx])); print(round(tab, 3))
}
cat("\nMSE relative to GEV MLE\n")
tab <- t(sapply(S, function(s) (s$mse / S[[1]]$mse)[idx])); rownames(tab) <- nm
colnames(tab) <- paste0("T=", round(Tp[idx])); print(round(tab, 3))
sink()
cat("\nwrote headline figures\n")

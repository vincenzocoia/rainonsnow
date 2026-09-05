# ---------------------------------------------------------------------------
# Summaries and figures from the Monte Carlo output.
#
# Outputs: out/fig-return-level.png, out/fig-mse.png, out/fig-decomposition.png,
#          out/fig-samplesize.png, out/results.txt
# ---------------------------------------------------------------------------
source("R/setup.R")
source("R/config.R")

METHODS <- c(mle = "GEV MLE", lmom = "GEV L-moments",
             cqe = "Composite quantile (L1)", cee = "Composite expectile (L2)")
COL <- c(mle = "#c0392b", lmom = "#e67e22", cqe = "#2980b9", cee = "#1e8449")
LTY <- c(mle = 2, lmom = 4, cqe = 3, cee = 1)

Tp <- RETURN_PERIODS
p_of_T <- 1 - 1 / Tp
truth_rl <- q_true(p_of_T)

# Return levels implied by each fitted parameter vector: N_REP x length(Tp).
return_levels <- function(par_mat) {
  t(apply(par_mat, 1, function(th)
    if (any(is.na(th))) rep(NA_real_, length(p_of_T))
    else qgev(p_of_T, th[1], th[2], th[3])))
}

summarise <- function(rl) {
  ok <- complete.cases(rl)
  rl <- rl[ok, , drop = FALSE]
  err <- sweep(rl, 2, truth_rl, "-")
  list(n_used = nrow(rl),
       mean = colMeans(rl), median = apply(rl, 2, median),
       bias = colMeans(err), var = apply(rl, 2, var),
       mse = colMeans(err^2),
       # Monte Carlo standard error of the MSE. Worth carrying: the sampling
       # distribution of a far-tail return level is heavy-tailed, so the MSE is
       # itself estimated with meaningful noise.
       mse_se = apply(err^2, 2, sd) / sqrt(nrow(err)),
       medae = apply(abs(err), 2, median),
       lo = apply(rl, 2, quantile, 0.025), hi = apply(rl, 2, quantile, 0.975))
}

load_n <- function(n) {
  d <- readRDS(sprintf("out/fits-n%d.rds", n))
  lapply(d$fits, function(f) summarise(return_levels(f)))
}

S <- list(); for (n in c(50, 100, 200)) S[[as.character(n)]] <- load_n(n)
s <- S[[as.character(N_OBS)]]

# --- figure 1: return period vs return level -------------------------------
png("out/fig-return-level.png", width = 1700, height = 1150, res = 135)
par(mfrow = c(2, 2), mar = c(4.3, 4.4, 3.2, 1.2))
ylim <- c(0, 45)
for (m in names(METHODS)) {
  plot(Tp, truth_rl, type = "n", log = "x", ylim = ylim,
       xlab = "return period T (years)", ylab = "return level",
       main = sprintf("%s  (n = %d)", METHODS[m], N_OBS))
  polygon(c(Tp, rev(Tp)), c(s[[m]]$lo, rev(s[[m]]$hi)),
          col = adjustcolor(COL[m], 0.20), border = NA)
  lines(Tp, s[[m]]$lo, col = COL[m], lwd = 1, lty = 3)
  lines(Tp, s[[m]]$hi, col = COL[m], lwd = 1, lty = 3)
  lines(Tp, s[[m]]$median, col = COL[m], lwd = 2.5)
  lines(Tp, truth_rl, col = "black", lwd = 2.5)
  legend("topleft", c("truth", "median estimate", "central 95% of estimates"),
         col = c("black", COL[m], adjustcolor(COL[m], 0.35)),
         lwd = c(2.5, 2.5, 8), bty = "n", cex = 0.82)
}
dev.off()

# --- figure 2: MSE over return period --------------------------------------
png("out/fig-mse.png", width = 1600, height = 700, res = 135)
par(mfrow = c(1, 2), mar = c(4.4, 4.6, 3.2, 1.2))
mse <- sapply(names(METHODS), function(m) s[[m]]$mse)
plot(Tp, mse[, 1], type = "n", log = "xy", ylim = range(mse, finite = TRUE),
     xlab = "return period T (years)", ylab = "MSE of estimated return level",
     main = sprintf("MSE over return period (n = %d)", N_OBS))
for (m in names(METHODS)) lines(Tp, s[[m]]$mse, col = COL[m], lwd = 2.6, lty = LTY[m])
legend("topleft", METHODS, col = COL, lwd = 2.6, lty = LTY, bty = "n", cex = 0.8)

rel <- sapply(names(METHODS), function(m) s[[m]]$mse / s$mle$mse)
plot(Tp, rel[, 1], type = "n", log = "x", ylim = c(0, max(2.2, min(4, max(rel)))),
     xlab = "return period T (years)", ylab = "MSE relative to GEV MLE",
     main = "MSE ratio (below 1 = better than MLE)")
abline(h = 1, col = "grey60")
for (m in names(METHODS)) lines(Tp, rel[, m], col = COL[m], lwd = 2.6, lty = LTY[m])
legend("topright", METHODS, col = COL, lwd = 2.6, lty = LTY, bty = "n", cex = 0.8)
dev.off()

# --- figure 3: bias^2 / variance decomposition ------------------------------
png("out/fig-decomposition.png", width = 1700, height = 620, res = 135)
par(mfrow = c(1, 3), mar = c(4.4, 4.6, 3.2, 1.2))
b2 <- sapply(names(METHODS), function(m) s[[m]]$bias^2)
vv <- sapply(names(METHODS), function(m) s[[m]]$var)
plot(Tp, b2[, 1], type = "n", log = "xy", ylim = range(b2[b2 > 0], finite = TRUE),
     xlab = "return period T", ylab = "bias^2", main = "Squared bias")
for (m in names(METHODS)) lines(Tp, s[[m]]$bias^2, col = COL[m], lwd = 2.6, lty = LTY[m])
legend("bottomright", METHODS, col = COL, lwd = 2.6, lty = LTY, bty = "n", cex = 0.75)
plot(Tp, vv[, 1], type = "n", log = "xy", ylim = range(vv, finite = TRUE),
     xlab = "return period T", ylab = "variance", main = "Variance")
for (m in names(METHODS)) lines(Tp, s[[m]]$var, col = COL[m], lwd = 2.6, lty = LTY[m])
plot(Tp, s$mle$medae, type = "n", log = "xy",
     ylim = range(sapply(names(METHODS), function(m) s[[m]]$medae)),
     xlab = "return period T", ylab = "median |error|",
     main = "Median absolute error\n(MSE-free view)")
for (m in names(METHODS)) lines(Tp, s[[m]]$medae, col = COL[m], lwd = 2.6, lty = LTY[m])
dev.off()

# --- figure 4: sample-size sensitivity --------------------------------------
png("out/fig-samplesize.png", width = 1700, height = 620, res = 135)
par(mfrow = c(1, 3), mar = c(4.4, 4.6, 3.2, 1.2))
for (n in c(50, 100, 200)) {
  ss <- S[[as.character(n)]]
  r <- sapply(names(METHODS), function(m) ss[[m]]$mse / ss$mle$mse)
  plot(Tp, r[, 1], type = "n", log = "x", ylim = c(0, 2.5),
       xlab = "return period T", ylab = "MSE relative to GEV MLE",
       main = sprintf("n = %d", n))
  abline(h = 1, col = "grey60")
  for (m in names(METHODS)) lines(Tp, r[, m], col = COL[m], lwd = 2.6, lty = LTY[m])
  if (n == 50) legend("topright", METHODS, col = COL, lwd = 2.6, lty = LTY,
                      bty = "n", cex = 0.72)
}
dev.off()

# --- numeric summary --------------------------------------------------------
sink("out/results.txt", split = TRUE)
cat("=== Monte Carlo results ===\n")
cat(sprintf("replicates: %d   truth: max(GEV(%g,%g,%g), N(%g,%g))\n\n",
            N_REP, TRUTH$mu, TRUTH$sigma, TRUTH$xi, TRUTH$a, TRUTH$b))
Tshow <- c(2, 5, 10, 20, 50, 100, 200, 500, 1000)
idx <- sapply(Tshow, function(t) which.min(abs(Tp - t)))
for (n in c(50, 100, 200)) {
  ss <- S[[as.character(n)]]
  cat(sprintf("\n---------------- n = %d ----------------\n", n))
  for (stat in c("bias", "sd", "rmse", "rel_rmse_pct", "mse_ratio_vs_mle",
                 "mse_MC_se_pct_of_mse", "medae")) {
    cat(sprintf("\n%s\n", stat))
    tab <- t(sapply(names(METHODS), function(m) {
      v <- switch(stat,
        bias = ss[[m]]$bias, sd = sqrt(ss[[m]]$var), rmse = sqrt(ss[[m]]$mse),
        rel_rmse_pct = 100 * sqrt(ss[[m]]$mse) / truth_rl,
        mse_ratio_vs_mle = ss[[m]]$mse / ss$mle$mse,
        mse_MC_se_pct_of_mse = 100 * ss[[m]]$mse_se / ss[[m]]$mse,
        medae = ss[[m]]$medae)
      v[idx]
    }))
    colnames(tab) <- paste0("T=", round(Tp[idx]))
    rownames(tab) <- METHODS
    print(round(tab, 3))
  }
}
sink()
cat("\nwrote figures and out/results.txt\n")

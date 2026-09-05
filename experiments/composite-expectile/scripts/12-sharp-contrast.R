# ---------------------------------------------------------------------------
# A setting where the composite expectile estimator wins outright.
#
# Same construction as the main study -- max(GEV, Normal), a product of cdfs --
# but with a much tighter, higher normal: N(2.5, 0.3) against GEV(0, 1, 0.2).
# Hydrologically this is the sharper version of the same story: in most years
# the annual maximum is an ordinary snowmelt peak, and those cluster tightly
# because they are set by a fairly repeatable snowpack; in the remaining years a
# rain-on-snow event produces something from a genuinely heavy-tailed process.
#
# Why it favours the composite estimators, where the earlier setting did not:
# a near-degenerate body drags a global fit hard (maximum likelihood ends up 44%
# low at the 1000-year level) while leaving the levels above the median entirely
# GEV-governed, so a tail-weighted criterion gets a clean and precise read. The
# fitted shape has a standard deviation of about 0.10 against a shape bias of
# roughly 0.30 in the global fits -- and that ratio, not the size of the bias
# alone, is what decides the contest.
#
# Output: out/fits-sharp.rds, out/sharp.txt, out/fig-sharp.png
# ---------------------------------------------------------------------------
source("R/setup.R")
source("R/config.R")
library(parallel)

TR <- list(mu = 0, sigma = 1, xi = 0.2, a = 2.5, b = 0.3)
P0_SET <- c(0.5, 0.7, 0.9)
n <- N_OBS
Tp <- RETURN_PERIODS; truth_rl <- q_true(1 - 1 / Tp, TR)
NC <- 2

# --- asymptotic targets -----------------------------------------------------
set.seed(606); big <- r_true(2e6, TR)
g50 <- make_level_grid(0.5, N_PANEL, N_GL); w50 <- make_weight(0.5, 0.53)
g90 <- make_level_grid(0.9, N_PANEL, N_GL); w90 <- make_weight(0.9, 0.93)
pseudo <- rbind(MLE = fit_mle(big), Lmoments = fit_lmom(big),
                `CQE p0=0.9` = fit_cqe(big, g90, w90),
                `CEE p0=0.5` = fit_cee(big, g50, w50))
colnames(pseudo) <- c("mu", "sigma", "xi")

# --- simulation -------------------------------------------------------------
set.seed(2718)
datasets <- lapply(seq_len(N_REP), function(i) r_true(n, TR))
cat("baselines ... "); t0 <- Sys.time()
b <- mclapply(datasets, function(y) { l <- fit_lmom(y); list(lmom = l, mle = fit_mle(y, start = l)) },
              mc.cores = NC)
fits <- list(mle = do.call(rbind, lapply(b, `[[`, "mle")),
             lmom = do.call(rbind, lapply(b, `[[`, "lmom")))
cat(sprintf("%.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
for (p0 in P0_SET) {
  g <- make_level_grid(p0, N_PANEL, N_GL); wf <- make_weight(p0, min(p0 + 0.03, 0.999))
  for (ty in c("quantile", "expectile")) {
    cat(sprintf("p0 = %.2f  %-9s ... ", p0, ty)); t0 <- Sys.time()
    fits[[sprintf("%s|%.2f", ty, p0)]] <- do.call(rbind, mclapply(datasets,
      function(y) fit_composite(y, g, wf, ty, start = fit_lmom(y)), mc.cores = NC))
    cat(sprintf("%.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  }
}
saveRDS(list(truth = TR, fits = fits, n = n, pseudo = pseudo), "out/fits-sharp.rds")

rls <- function(m) t(apply(m, 1, function(th) if (any(is.na(th))) rep(NA_real_, length(Tp))
                           else qgev(1 - 1 / Tp, th[1], th[2], th[3])))
R <- lapply(fits, rls)
summ <- function(r) { r <- r[complete.cases(r), , drop = FALSE]
  e <- sweep(r, 2, truth_rl, "-")
  list(bias = colMeans(e), sd = apply(r, 2, sd), mse = colMeans(e^2)) }
S <- lapply(R, summ)
Tshow <- c(10, 20, 50, 100, 200, 500, 1000)
idx <- sapply(Tshow, function(t) which.min(abs(Tp - t)))

sink("out/sharp.txt", split = TRUE)
cat(sprintf("=== Sharp contrast: max(GEV(%g,%g,%g), N(%g,%g)) ===\n",
            TR$mu, TR$sigma, TR$xi, TR$a, TR$b))
cat(sprintf("n = %d, %d replicates\n\n", n, N_REP))
cat("--- asymptotic targets (fit to n = 2e6); true tail (0, 1, 0.2) ---\n")
print(round(pseudo, 4))
cat("\nasymptotic % bias in the return level\n")
ab <- t(apply(pseudo, 1, function(th)
  100 * (qgev(1 - 1 / Tp[idx], th[1], th[2], th[3]) - truth_rl[idx]) / truth_rl[idx]))
colnames(ab) <- paste0("T=", round(Tp[idx])); print(round(ab, 1))
cat("\n--- median fitted shape and its sd (true xi = 0.2) ---\n")
print(round(data.frame(median = sapply(fits, function(m) median(m[, 3], na.rm = TRUE)),
                       sd = sapply(fits, function(m) sd(m[, 3], na.rm = TRUE))), 3))
for (lab in c("mse ratio vs GEV MLE", "bias", "sd")) {
  cat(sprintf("\n--- %s ---\n", lab))
  tab <- t(sapply(S, function(s) switch(lab, `mse ratio vs GEV MLE` = s$mse / S$mle$mse,
                                        bias = s$bias, sd = s$sd)[idx]))
  colnames(tab) <- paste0("T=", round(Tp[idx])); print(round(tab, 3))
}
cat("\n--- paired MC standard error of each MSE ratio ---\n")
se <- t(sapply(names(R), function(nm) {
  ok <- complete.cases(R[[nm]]) & complete.cases(R$mle)
  ea <- sweep(R[[nm]][ok, ], 2, truth_rl, "-")^2
  em <- sweep(R$mle[ok, ], 2, truth_rl, "-")^2
  (apply(ea - em, 2, sd) / sqrt(sum(ok)) / colMeans(em))[idx] }))
colnames(se) <- paste0("T=", round(Tp[idx])); print(round(se, 3))
sink()

png("out/fig-sharp.png", width = 1650, height = 640, res = 133)
par(mfrow = c(1, 3), mar = c(4.4, 4.6, 3.2, 1.2))
xs <- seq(-1, 12, length.out = 1200)
plot(xs, d_true(xs, TR), type = "l", lwd = 2.8, xlab = "x", ylab = "density",
     main = "Truth: a tight body, a heavy tail")
lines(xs, dgev(xs, TR$mu, TR$sigma, TR$xi), lwd = 2, col = "#c0392b", lty = 2)
lines(xs, dnorm(xs, TR$a, TR$b), lwd = 1.8, col = "#2980b9", lty = 3)
legend("topright", c("truth = max(GEV, N)", "GEV component", "normal component"),
       col = c("black", "#c0392b", "#2980b9"), lty = c(1, 2, 3), lwd = 2, bty = "n", cex = 0.75)
cols <- colorRampPalette(c("#a9cce3", "#154360"))(length(P0_SET))
plot(Tp, rep(1, length(Tp)), type = "n", log = "xy", ylim = c(0.3, 30),
     xlab = "return period T (years)", ylab = "MSE relative to GEV MLE",
     main = "MSE relative to GEV MLE")
abline(h = 1, col = "grey50", lwd = 1.5)
lines(Tp, S$lmom$mse / S$mle$mse, col = "#e67e22", lwd = 3, lty = 2)
for (i in seq_along(P0_SET)) {
  lines(Tp, S[[sprintf("expectile|%.2f", P0_SET[i])]]$mse / S$mle$mse, col = "#1e8449",
        lwd = 2.6, lty = i)
  lines(Tp, S[[sprintf("quantile|%.2f", P0_SET[i])]]$mse / S$mle$mse, col = "#2980b9",
        lwd = 2.6, lty = i)
}
legend("topleft", c("L2 (expectile)", "L1 (quantile)", "GEV L-moments",
                    sprintf("p0 = %.1f", P0_SET)),
       col = c("#1e8449", "#2980b9", "#e67e22", rep("grey40", 3)),
       lwd = c(2.6, 2.6, 3, 1, 1, 1), lty = c(1, 1, 2, 1, 2, 3), bty = "n", cex = 0.72)
plot(Tp, truth_rl, type = "n", log = "x", ylim = c(0, max(truth_rl) * 1.5),
     xlab = "return period T (years)", ylab = "return level", main = "Median fitted return level")
lines(Tp, truth_rl, col = "black", lwd = 3)
for (nm in c("mle", "lmom", "expectile|0.50", "quantile|0.90"))
  lines(Tp, apply(R[[nm]], 2, median, na.rm = TRUE), lwd = 2.4,
        col = switch(nm, mle = "#c0392b", lmom = "#e67e22",
                     `expectile|0.50` = "#1e8449", `quantile|0.90` = "#2980b9"))
legend("topleft", c("truth", "GEV MLE", "GEV L-moments", "L2, p0 = 0.5", "L1, p0 = 0.9"),
       col = c("black", "#c0392b", "#e67e22", "#1e8449", "#2980b9"), lwd = 2.4, bty = "n", cex = 0.75)
dev.off()
cat("\nwrote out/fig-sharp.png and out/sharp.txt\n")

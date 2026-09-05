# ---------------------------------------------------------------------------
# A setup where the global fits MUST fail: a grafted truth.
#
# max(GEV, Normal) turned out to be a gentle test. Because the GEV component
# still shapes the upper body, the body leaks tail information into a global
# fit, and maximum likelihood recovers enough of the tail to stay competitive.
#
# A graft removes that leak. Above a join level the truth IS the GEV, exactly;
# below it the body is a different distribution carrying no information about
# the tail at all. Hydrologically this is the standard mixed-population picture,
# and it fits the rain-on-snow story directly: in most years the annual maximum
# is an ordinary snowmelt event, drawn from a much less variable and physically
# limited distribution; in the remaining years it is a rain-on-snow event, and
# those are what the GEV tail describes.
#
#   truth: GEV(0, 1, 0.4) above the 90th percentile of that GEV (u = 3.65),
#          N(0.65, 1.58) truncated above u below it, density matched at the join.
#
# Asymptotically this leaves maximum likelihood 55% low at the 1000-year level
# and L-moments 45% low, against under 2% for both composite estimators.
#
# Output: out/fits-graft.rds, out/graft.txt, out/fig-graft-dgp.png, out/fig-graft.png
# ---------------------------------------------------------------------------
source("R/setup.R")
source("R/config.R")
library(parallel)

TR <- make_graft(mu = 0, sigma = 1, xi = 0.4, p_join = 0.90, a_offset = 3, branch = "tight")
P0_SET <- c(0.5, 0.7, 0.9)
n <- N_OBS
Tp <- RETURN_PERIODS; truth_rl <- q_graft(1 - 1 / Tp, TR)

# --- picture of the truth ---------------------------------------------------
png("out/fig-graft-dgp.png", width = 1600, height = 620, res = 133)
par(mfrow = c(1, 2), mar = c(4.3, 4.4, 3.2, 1.2))
xs <- seq(TR$a - 4 * TR$b, q_graft(0.995, TR), length.out = 1400)
plot(xs, d_graft(xs, TR), type = "l", lwd = 2.8, col = "black",
     xlab = "x", ylab = "density", main = "Grafted truth: density")
lines(xs, dgev(xs, TR$mu, TR$sigma, TR$xi), lwd = 2, col = "#c0392b", lty = 2)
abline(v = TR$u, col = "#27ae60", lty = 3, lwd = 2)
legend("topright", c("truth", "GEV tail component", "join"),
       col = c("black", "#c0392b", "#27ae60"), lty = c(1, 2, 3), lwd = 2, bty = "n", cex = 0.8)
Tg <- exp(seq(log(1.05), log(2000), length.out = 400)); pg <- 1 - 1 / Tg
plot(Tg, q_graft(pg, TR), type = "l", lwd = 2.8, log = "x",
     xlab = "return period T (years)", ylab = "return level",
     main = "Grafted truth: frequency-magnitude")
lines(Tg, qgev(pg, TR$mu, TR$sigma, TR$xi), lwd = 2, col = "#c0392b", lty = 2)
abline(v = 1 / (1 - TR$Fu), col = "#27ae60", lty = 3, lwd = 2)
legend("topleft", c("truth", "GEV tail component", "join (T = 10)"),
       col = c("black", "#c0392b", "#27ae60"), lty = c(1, 2, 3), lwd = 2, bty = "n", cex = 0.8)
dev.off()

# --- simulation -------------------------------------------------------------
set.seed(31415)
datasets <- lapply(seq_len(N_REP), function(i) r_graft(n, TR))
cat("baselines ... "); t0 <- Sys.time()
b <- mclapply(datasets, function(y) { l <- fit_lmom(y); list(lmom = l, mle = fit_mle(y, start = l)) },
              mc.cores = 2)
fits <- list(mle = do.call(rbind, lapply(b, `[[`, "mle")),
             lmom = do.call(rbind, lapply(b, `[[`, "lmom")))
cat(sprintf("%.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
for (p0 in P0_SET) {
  g <- make_level_grid(p0, N_PANEL, N_GL); wf <- make_weight(p0, min(p0 + 0.03, 0.999))
  for (ty in c("quantile", "expectile")) {
    cat(sprintf("p0 = %.2f  %-9s ... ", p0, ty)); t0 <- Sys.time()
    fits[[sprintf("%s|%.2f", ty, p0)]] <- do.call(rbind, mclapply(datasets,
      function(y) fit_composite(y, g, wf, ty, start = fit_lmom(y)), mc.cores = 2))
    cat(sprintf("%.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  }
}
saveRDS(list(truth = TR, fits = fits, n = n), "out/fits-graft.rds")

rls <- function(m) t(apply(m, 1, function(th) if (any(is.na(th))) rep(NA_real_, length(Tp))
                           else qgev(1 - 1 / Tp, th[1], th[2], th[3])))
R <- lapply(fits, rls)
summ <- function(r) { r <- r[complete.cases(r), , drop = FALSE]
  e <- sweep(r, 2, truth_rl, "-")
  list(bias = colMeans(e), sd = apply(r, 2, sd), mse = colMeans(e^2)) }
S <- lapply(R, summ)
Tshow <- c(10, 20, 50, 100, 200, 500, 1000)
idx <- sapply(Tshow, function(t) which.min(abs(Tp - t)))

sink("out/graft.txt", split = TRUE)
cat(sprintf("=== Grafted truth: GEV(%g, %g, %g) above u = %.3f, N(%.3f, %.3f) truncated below ===\n",
            TR$mu, TR$sigma, TR$xi, TR$u, TR$a, TR$b))
cat(sprintf("join at the %.0fth percentile; n = %d, %d replicates\n\n", 100 * TR$Fu, n, N_REP))
cat("--- median fitted shape (true tail xi = 0.4) ---\n")
print(round(sapply(fits, function(m) median(m[, 3], na.rm = TRUE)), 3))
for (lab in c("mse ratio vs GEV MLE", "bias", "sd")) {
  cat(sprintf("\n--- %s ---\n", lab))
  tab <- t(sapply(S, function(s) switch(lab, `mse ratio vs GEV MLE` = s$mse / S$mle$mse,
                                        bias = s$bias, sd = s$sd)[idx]))
  colnames(tab) <- paste0("T=", round(Tp[idx])); print(round(tab, 3))
}
cat("\n--- paired MC standard error of the MSE difference from the MLE, as a ratio ---\n")
se <- t(sapply(names(R), function(nm) {
  ok <- complete.cases(R[[nm]]) & complete.cases(R$mle)
  ea <- sweep(R[[nm]][ok, ], 2, truth_rl, "-")^2
  em <- sweep(R$mle[ok, ], 2, truth_rl, "-")^2
  (apply(ea - em, 2, sd) / sqrt(sum(ok)) / colMeans(em))[idx] }))
colnames(se) <- paste0("T=", round(Tp[idx])); print(round(se, 3))
sink()

png("out/fig-graft.png", width = 1650, height = 640, res = 133)
par(mfrow = c(1, 3), mar = c(4.4, 4.6, 3.2, 1.2))
cols <- colorRampPalette(c("#a9cce3", "#154360"))(length(P0_SET))
for (ty in c("quantile", "expectile")) {
  plot(Tp, rep(1, length(Tp)), type = "n", log = "xy", ylim = c(0.05, 30),
       xlab = "return period T (years)", ylab = "MSE relative to GEV MLE",
       main = sprintf("Composite %s (%s)", ty, if (ty == "quantile") "L1" else "L2"))
  abline(h = 1, col = "grey50", lwd = 1.5)
  lines(Tp, S$lmom$mse / S$mle$mse, col = "#e67e22", lwd = 3, lty = 2)
  for (i in seq_along(P0_SET))
    lines(Tp, S[[sprintf("%s|%.2f", ty, P0_SET[i])]]$mse / S$mle$mse, col = cols[i], lwd = 2.6)
  legend("bottomleft", c(sprintf("p0 = %.2f", P0_SET), "GEV L-moments"),
         col = c(cols, "#e67e22"), lwd = 2.6, lty = c(rep(1, length(P0_SET)), 2),
         bty = "n", cex = 0.72)
}
plot(Tp, truth_rl, type = "n", log = "x", ylim = c(0, max(truth_rl) * 1.6),
     xlab = "return period T (years)", ylab = "return level",
     main = "Median fitted return level")
lines(Tp, truth_rl, col = "black", lwd = 3)
for (nm in c("mle", "lmom", "expectile|0.50", "quantile|0.90")) {
  m <- apply(R[[nm]], 2, median, na.rm = TRUE)
  lines(Tp, m, lwd = 2.4, col = switch(nm, mle = "#c0392b", lmom = "#e67e22",
                                       `expectile|0.50` = "#1e8449", `quantile|0.90` = "#2980b9"))
}
legend("topleft", c("truth", "GEV MLE", "GEV L-moments", "L2, p0 = 0.5", "L1, p0 = 0.9"),
       col = c("black", "#c0392b", "#e67e22", "#1e8449", "#2980b9"), lwd = 2.4, bty = "n", cex = 0.75)
dev.off()
cat("\nwrote out/fig-graft.png, out/fig-graft-dgp.png, out/graft.txt\n")

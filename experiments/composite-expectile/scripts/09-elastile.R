# ---------------------------------------------------------------------------
# The alpha-elastile estimator: does mixing the two losses beat either alone?
#
# Part 1 (population): how much of the body contamination survives into the
#   alpha-elastile, and what each alpha converges to. This predicts the bias.
# Part 2 (simulation): the alpha sweep at the primary sample size.
#
# Output: out/fits-elastile.rds, out/elastile.txt, out/fig-elastile.png
# ---------------------------------------------------------------------------
source("R/setup.R")
source("R/config.R")
library(parallel)

ALPHAS <- c(0, 0.02, 0.05, 0.1, 0.2, 0.35, 0.5, 0.7, 0.85, 1)
tr <- TRUTH

# --- the scale s at the population level, same convention as the estimator ---
risk_L1 <- function(p, q, xi) gev_partial_moment(q, 0, 1, xi) + (1 - p) * (q - gev_mean(0, 1, xi))
risk_L2 <- function(p, e, xi) {
  m <- gev_mean(0, 1, xi); v <- (gamma(1 - 2 * xi) - gamma(1 - xi)^2) / xi^2
  up <- gev_second_partial_moment(e, 0, 1, xi)
  p * up + (1 - p) * (v + (m - e)^2 - up)
}
pg <- GRID$p; wg <- GRID$w_quad * W_FUN(pg)
S_POP <- sum(wg * risk_L2(pg, gev_expectile(pg, 0, 1, tr$xi), tr$xi)) /
         sum(wg * risk_L1(pg, qgev(pg, 0, 1, tr$xi), tr$xi))

sink("out/elastile.txt", split = TRUE)
cat("=== The alpha-elastile ===\n")
cat("At level p it is the root of\n\n")
cat("  G(t) = (2 alpha / s) [ (1 - 2p) phi(t) + (1 - p)(t - m) ] + (1 - alpha) [ F(t) - p ]\n\n")
cat("The first bracket vanishes at the expectile, the second at the quantile, and\n")
cat("both are increasing in t, so G is strictly increasing: the root is unique and\n")
cat("always lies between the quantile and the expectile. alpha = 0 gives the\n")
cat("quantile, alpha = 1 the expectile.\n\n")
cat(sprintf("scale s (population convention) = %.4f\n\n", S_POP))

cat("=== 1. How much body contamination survives into the alpha-elastile ===\n")
cat("Percent by which the TRUTH's functional exceeds its GEV component's.\n")
cat("At alpha = 0 this is the quantile contamination (dies out fast); at alpha = 1\n")
cat("the expectile contamination (does not).\n\n")
pp <- c(0.95, 0.98, 0.99, 0.999)
m_tr <- mean_true(tr)
tab <- t(sapply(ALPHAS, function(a) {
  tt <- elastile_true(pp, a, S_POP, tr, m_tr)
  gg <- if (a <= 0) qgev(pp, tr$mu, tr$sigma, tr$xi)
        else if (a >= 1) gev_expectile(pp, tr$mu, tr$sigma, tr$xi)
        else gev_elastile(pp, tr$mu, tr$sigma, tr$xi, a, S_POP)
  100 * (tt - gg) / gg
}))
rownames(tab) <- sprintf("alpha=%.2f", ALPHAS)
colnames(tab) <- sprintf("p=%.3f (T=%.0f)", pp, 1 / (1 - pp))
print(round(tab, 3))

cat("\n=== 2. Asymptotic target of each alpha (fit to n = 2e6) ===\n")
set.seed(20260905); big <- r_true(2e6)
Tp <- RETURN_PERIODS; truth_rl <- q_true(1 - 1 / Tp)
Tshow <- c(10, 20, 50, 100, 200, 500, 1000)
idx <- sapply(Tshow, function(t) which.min(abs(Tp - t)))
pseudo <- t(sapply(ALPHAS, function(a) fit_elastile(big, GRID, W_FUN, a)))
rownames(pseudo) <- sprintf("alpha=%.2f", ALPHAS)
colnames(pseudo) <- c("mu", "sigma", "xi")
print(round(pseudo, 4))
cat(sprintf("\ntrue GEV tail: mu = %g, sigma = %g, xi = %g\n", tr$mu, tr$sigma, tr$xi))
cat("\nasymptotic % bias in the return level\n")
ab <- t(apply(pseudo, 1, function(th) 100 * (qgev(1 - 1 / Tp[idx], th[1], th[2], th[3]) -
                                             truth_rl[idx]) / truth_rl[idx]))
colnames(ab) <- paste0("T=", round(Tp[idx])); print(round(ab, 2))
sink()

# --- Part 2: the simulation sweep -------------------------------------------
#
# Run at two weights: the one the experiment was designed around (p0 = 0.95,
# placed where the QUANTILE contamination dies out) and the one the tail-focus
# sweep found to be much better for the expectile side (p0 = 0.5). Asking
# whether mixing helps is only meaningful where each pure loss is at its best.
n <- N_OBS
WEIGHTS <- list(`p0=0.50` = 0.50, `p0=0.95` = 0.95)
set.seed(4242 + n)                              # the same datasets as scripts/01
datasets <- lapply(seq_len(N_REP), function(i) r_true(n))
fits <- list()
for (wn in names(WEIGHTS)) {
  p0 <- WEIGHTS[[wn]]
  g <- make_level_grid(p0, N_PANEL, N_GL); wf <- make_weight(p0, min(p0 + 0.03, 0.999))
  for (a in ALPHAS) {
    cat(sprintf("%s alpha = %.2f ... ", wn, a)); t0 <- Sys.time()
    fits[[sprintf("%s|%.2f", wn, a)]] <- do.call(rbind, mclapply(datasets,
      function(y) fit_elastile(y, g, wf, a, start = fit_lmom(y)),
      mc.cores = detectCores()))
    cat(sprintf("%.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  }
}
saveRDS(list(alphas = ALPHAS, weights = WEIGHTS, fits = fits, n = n, s_pop = S_POP),
        "out/fits-elastile.rds")

summ <- function(m) {
  rl <- t(apply(m, 1, function(th) if (any(is.na(th))) rep(NA_real_, length(Tp))
                else qgev(1 - 1 / Tp, th[1], th[2], th[3])))
  rl <- rl[complete.cases(rl), , drop = FALSE]
  e <- sweep(rl, 2, truth_rl, "-")
  list(bias = colMeans(e), sd = apply(rl, 2, sd), mse = colMeans(e^2),
       xi = median(m[, 3], na.rm = TRUE))
}
S <- lapply(fits, summ)
base <- readRDS(sprintf("out/fits-n%d.rds", n))$fits
mle <- summ(base$mle); lmom <- summ(base$lmom)

sink("out/elastile.txt", append = TRUE, split = TRUE)
cat(sprintf("\n\n=== 3. Simulation, n = %d, %d replicates ===\n", n, N_REP))
for (wn in names(WEIGHTS)) {
  Sw <- S[grep(paste0("^", wn), names(S))]
  cat(sprintf("\n################  %s  ################\n", wn))
  for (lab in c("mse ratio vs GEV MLE", "bias", "sd")) {
    cat(sprintf("\n--- %s ---\n", lab))
    tab <- t(sapply(Sw, function(s) switch(lab, `mse ratio vs GEV MLE` = s$mse / mle$mse,
                                           bias = s$bias, sd = s$sd)[idx]))
    rownames(tab) <- sprintf("alpha=%.2f", ALPHAS)
    ref <- rbind(`GEV MLE` = switch(lab, `mse ratio vs GEV MLE` = rep(1, length(idx)),
                                    bias = mle$bias[idx], sd = mle$sd[idx]),
                 `GEV L-moments` = switch(lab, `mse ratio vs GEV MLE` = (lmom$mse / mle$mse)[idx],
                                          bias = lmom$bias[idx], sd = lmom$sd[idx]))
    tab <- rbind(tab, ref); colnames(tab) <- paste0("T=", round(Tp[idx]))
    print(round(tab, 3))
  }
  cat("\n--- median fitted shape (true xi = 0.2) ---\n")
  print(round(setNames(sapply(Sw, `[[`, "xi"), sprintf("a=%.2f", ALPHAS)), 3))
  cat("\n--- alpha with the lowest MSE at each return period ---\n")
  best <- sapply(idx, function(k) ALPHAS[which.min(sapply(Sw, function(s) s$mse[k]))])
  print(setNames(best, paste0("T=", round(Tp[idx]))))
  cat("--- and how it compares with the better of alpha = 0 and alpha = 1 ---\n")
  gain <- sapply(idx, function(k) {
    m <- sapply(Sw, function(s) s$mse[k]); min(m) / min(m[1], m[length(m)]) })
  print(round(setNames(gain, paste0("T=", round(Tp[idx]))), 3))
}
sink()

png("out/fig-elastile.png", width = 1650, height = 1080, res = 134)
par(mfrow = c(2, 3), mar = c(4.4, 4.6, 3.4, 1.2))
cols <- colorRampPalette(c("#2980b9", "#7d3c98", "#1e8449"))(length(ALPHAS))
for (wn in names(WEIGHTS)) {
  Sw <- S[grep(paste0("^", wn), names(S))]
  plot(Tp, rep(1, length(Tp)), type = "n", log = "xy", ylim = c(0.5, 40),
       xlab = "return period T", ylab = "MSE relative to GEV MLE",
       main = sprintf("MSE by alpha, %s (n = %d)", wn, n))
  abline(h = 1, col = "grey50", lwd = 1.5)
  lines(Tp, lmom$mse / mle$mse, col = "#e67e22", lwd = 3, lty = 2)
  for (i in seq_along(ALPHAS)) lines(Tp, Sw[[i]]$mse / mle$mse, col = cols[i], lwd = 2.3)
  legend("topright", c(sprintf("alpha = %.2f", ALPHAS), "L-moments"),
         col = c(cols, "#e67e22"), lwd = 2.3, lty = c(rep(1, length(ALPHAS)), 2),
         bty = "n", cex = 0.58)
  for (k in c(which.min(abs(Tp - 100)), which.min(abs(Tp - 1000)))) {
    r <- sapply(Sw, function(s) s$mse[k]) / mle$mse[k]
    plot(ALPHAS, r, type = "b", pch = 19, lwd = 2.4, col = "#7d3c98", log = "y",
         ylim = range(c(r, 1, lmom$mse[k] / mle$mse[k])),
         xlab = "alpha  (0 = quantile, 1 = expectile)",
         ylab = "MSE relative to GEV MLE",
         main = sprintf("%s,  T = %.0f", wn, Tp[k]))
    abline(h = 1, col = "#c0392b", lwd = 2)
    abline(h = lmom$mse[k] / mle$mse[k], col = "#e67e22", lwd = 2, lty = 2)
    legend("topright", c("elastile", "GEV MLE", "GEV L-moments"),
           col = c("#7d3c98", "#c0392b", "#e67e22"), lwd = 2.4, lty = c(1, 1, 2),
           bty = "n", cex = 0.72)
  }
}
dev.off()
cat("\nwrote out/fig-elastile.png and out/elastile.txt\n")

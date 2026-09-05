# ---------------------------------------------------------------------------
# How much of the result is the tail focus, and how much is quantiles vs
# expectiles?
#
# Sweeps the lower edge of the weight from p0 = 0 (w == 1 everywhere, no tail
# focus at all) up to p0 = 0.95. At p0 = 0 the composite quantile estimator is
# exactly CRPS minimisation, which is the natural benchmark against L-moments:
# both are global, quantile-flavoured, and neither has any tail focus.
#
# Output: out/fits-tailfocus.rds, out/tail-focus.txt, out/fig-tail-focus.png
# ---------------------------------------------------------------------------
source("R/setup.R")
source("R/config.R")
library(parallel)

P0_SET <- c(0, 0.5, 0.7, 0.8, 0.9, 0.95)
n <- N_OBS
set.seed(4242 + n)                              # the same datasets as scripts/01
datasets <- lapply(seq_len(N_REP), function(i) r_true(n))

wfun <- function(p0) {
  if (p0 <= 0) return(function(p) rep(1, length(p)))   # w == 1: no tail focus
  make_weight(p0, min(p0 + 0.03, 0.999))
}

fits <- list()
for (p0 in P0_SET) {
  g <- make_level_grid(p0, N_PANEL, N_GL); wf <- wfun(p0)
  for (ty in c("quantile", "expectile")) {
    cat(sprintf("p0 = %.2f  %-9s ... ", p0, ty)); t0 <- Sys.time()
    fits[[sprintf("%s|%.2f", ty, p0)]] <- do.call(rbind, mclapply(datasets,
      function(y) fit_composite(y, g, wf, ty, start = fit_lmom(y)),
      mc.cores = detectCores()))
    cat(sprintf("%.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  }
}
saveRDS(list(p0 = P0_SET, fits = fits, n = n), "out/fits-tailfocus.rds")

Tp <- RETURN_PERIODS; truth_rl <- q_true(1 - 1 / Tp)
summ <- function(m) {
  rl <- t(apply(m, 1, function(th) if (any(is.na(th))) rep(NA_real_, length(Tp))
                else qgev(1 - 1 / Tp, th[1], th[2], th[3])))
  rl <- rl[complete.cases(rl), , drop = FALSE]
  e <- sweep(rl, 2, truth_rl, "-")
  list(bias = colMeans(e), sd = apply(rl, 2, sd), mse = colMeans(e^2),
       xi = median(m[, 3], na.rm = TRUE),
       bound = mean(m[, 3] <= XI_LO + 1e-6, na.rm = TRUE))
}
S <- lapply(fits, summ)
base <- readRDS(sprintf("out/fits-n%d.rds", n))$fits
mle <- summ(base$mle); lmom <- summ(base$lmom)

sink("out/tail-focus.txt", split = TRUE)
cat(sprintf("=== Tail focus swept from none to strong, n = %d, %d replicates ===\n", n, N_REP))
cat("p0 = 0 means w == 1: no tail focus. For the quantile loss that is exactly\n")
cat("CRPS minimisation -- the like-for-like comparison against L-moments.\n\n")
Tshow <- c(10, 20, 50, 100, 200, 500, 1000)
idx <- sapply(Tshow, function(t) which.min(abs(Tp - t)))
cat("--- median fitted shape (true xi = 0.2), and fraction on the shape bound ---\n")
print(round(data.frame(xi = sapply(S, `[[`, "xi"), at_bound = sapply(S, `[[`, "bound")), 3))
cat(sprintf("\nfor reference: MLE xi = %.3f, L-moments xi = %.3f\n",
            median(base$mle[, 3], na.rm = TRUE), median(base$lmom[, 3], na.rm = TRUE)))
for (lab in c("mse ratio vs GEV MLE", "bias", "sd")) {
  cat(sprintf("\n--- %s ---\n", lab))
  tab <- t(sapply(S, function(s) switch(lab, `mse ratio vs GEV MLE` = s$mse / mle$mse,
                                        bias = s$bias, sd = s$sd)[idx]))
  ref <- rbind(`GEV MLE` = switch(lab, `mse ratio vs GEV MLE` = rep(1, length(idx)),
                                  bias = mle$bias[idx], sd = mle$sd[idx]),
               `GEV L-moments` = switch(lab, `mse ratio vs GEV MLE` = (lmom$mse / mle$mse)[idx],
                                        bias = lmom$bias[idx], sd = lmom$sd[idx]))
  tab <- rbind(tab, ref); colnames(tab) <- paste0("T=", round(Tp[idx]))
  print(round(tab, 3))
}
sink()

png("out/fig-tail-focus.png", width = 1600, height = 660, res = 134)
par(mfrow = c(1, 2), mar = c(4.4, 4.6, 3.2, 1.2))
cols <- colorRampPalette(c("#d5d8dc", "#154360"))(length(P0_SET))
for (ty in c("quantile", "expectile")) {
  plot(Tp, rep(1, length(Tp)), type = "n", log = "xy", ylim = c(0.3, 60),
       xlab = "return period T (years)", ylab = "MSE relative to GEV MLE",
       main = sprintf("Composite %s (%s), n = %d", ty,
                      if (ty == "quantile") "L1" else "L2", n))
  abline(h = 1, col = "grey50", lwd = 1.5)
  lines(Tp, lmom$mse / mle$mse, col = "#e67e22", lwd = 3, lty = 2)
  for (i in seq_along(P0_SET))
    lines(Tp, S[[sprintf("%s|%.2f", ty, P0_SET[i])]]$mse / mle$mse, col = cols[i], lwd = 2.6)
  legend("topright", c(sprintf("p0 = %.2f", P0_SET), "GEV L-moments"),
         col = c(cols, "#e67e22"), lwd = c(rep(2.6, length(P0_SET)), 3),
         lty = c(rep(1, length(P0_SET)), 2), bty = "n", cex = 0.72)
}
dev.off()
cat("\nwrote out/fig-tail-focus.png and out/tail-focus.txt\n")

# ---------------------------------------------------------------------------
# A self-contained GPD study.
#
# An astute practitioner facing a body that is obviously not GEV would not fit a
# single parametric distribution at all. They would pick a threshold, fit a GPD
# to the exceedances by maximum likelihood or L-moments, and hard graft it onto
# an empirical body. That is the competitor here, at three thresholds, and it is
# the reference against which everything is measured -- no GEV method appears in
# this comparison.
#
# The composite alternative fits a three-parameter GPD to the whole sample with
# a weight that rises from p = 0, then smooth grafts it onto the empirical body
# using that same weight as the handover. Unlike the GEV case, tying the two
# weights together is natural here: the composite weight already starts at zero
# on the body, which is exactly where the empirical distribution should be
# trusted, so the two jobs coincide.
#
# One caveat on the truth: max(GEV, Normal) has a tail that is exactly GEV, and
# a GEV tail is only ASYMPTOTICALLY GPD (S = 1 - exp(-z) against z, a relative
# error of about z/2, so ~5% at the 90th percentile and ~0.5% at the 99th). The
# GPD is therefore mildly misspecified for every method in this comparison,
# which is realistic and affects them all alike.
#
# Output: out/fits-gpd.rds, out/gpd.txt, out/fig-gpd.png
# ---------------------------------------------------------------------------
.libPaths(c("/home/user/Rlib-graft", .libPaths()))
source("R/setup.R"); source("R/config.R")
source("R/gpd.R"); source("R/gpd_estimators.R"); source("R/smoothgraft.R")
suppressMessages({library(distionary); library(distplyr)})
library(parallel)

n <- N_OBS; NC <- detectCores()
Tp <- RETURN_PERIODS; ex <- 1 / Tp; truth_rl <- q_true(1 - ex)
set.seed(4242 + n)
datasets <- lapply(seq_len(N_REP), function(i) r_true(n))

RL <- list()
# --- the practitioner's competitor: POT + hard graft -------------------------
for (v in c(0.85, 0.90, 0.95)) {
  cat(sprintf("POT-MLE  threshold %.2f ... ", v)); t0 <- Sys.time()
  RL[[sprintf("POT-MLE u=%.2f", v)]] <- do.call(rbind, mclapply(datasets,
    function(y) pot_return_levels(y, fit_pot_mle(y, v), ex), mc.cores = NC))
  cat(sprintf("%.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
cat("POT-Lmom threshold 0.90 ... ")
RL[["POT-Lmom u=0.90"]] <- do.call(rbind, mclapply(datasets,
  function(y) pot_return_levels(y, fit_pot_lmom(y, 0.90), ex), mc.cores = NC))
cat("done\n")
RL[["empirical"]] <- do.call(rbind, mclapply(datasets, function(y)
  as.numeric(eval_quantile(dst_empirical(y), at = 1 - ex)), mc.cores = NC))

# --- the composite alternative: composite GPD + smooth graft -----------------
# Weights that all rise from p = 0, differing in how fast. The smoothstep is
# already substantially on by mid-distribution; the powers rise immediately but
# stay small on the body, which is closer to "starts rising right away" while
# still concentrating on the tail.
COMP <- list(
  `L2, w = smoothstep(p/0.8)` = list(type = "expectile", w = make_weight(0, 0.80), anch = FALSE),
  `L2, w = p^2`               = list(type = "expectile", w = function(p) p^2,      anch = FALSE),
  `L2, w = p^6`               = list(type = "expectile", w = function(p) p^6,      anch = FALSE),
  `L1, w = p^2`               = list(type = "quantile",  w = function(p) p^2,      anch = FALSE)
)
grid0 <- make_level_grid(0, N_PANEL, N_GL)
PAR <- list()
for (nm in names(COMP)) {
  cc <- COMP[[nm]]
  cat(sprintf("fit %-28s ... ", nm)); t0 <- Sys.time()
  PAR[[nm]] <- do.call(rbind, mclapply(datasets,
    function(y) fit_gpd_composite(y, grid0, cc$w, cc$type, cc$anch), mc.cores = NC))
  cat(sprintf("%.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  RL[[paste0(nm, " (ungrafted)")]] <- t(apply(PAR[[nm]], 1, function(th)
    if (any(is.na(th))) rep(NA_real_, length(Tp)) else qgpd(1 - ex, th[1], th[2], th[3])))
  cat(sprintf("graft %-26s ... ", nm)); t0 <- Sys.time()
  RL[[paste0(nm, " + graft")]] <- do.call(rbind, mclapply(seq_along(datasets),
    function(i) as.numeric(graft_return_levels(datasets[[i]], PAR[[nm]][i, ],
                                               NA, NA, ex, family = "gpd", w_fun = cc$w)),
    mc.cores = NC))
  cat(sprintf("%.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}

# --- asymptotic targets ------------------------------------------------------
set.seed(101); big <- r_true(3e5)
asym <- list()
for (v in c(0.85, 0.90, 0.95)) {
  f <- fit_pot_mle(big, v)
  asym[[sprintf("POT-MLE u=%.2f", v)]] <-
    100 * (pot_return_levels(big, f, ex) - truth_rl) / truth_rl
}
f <- fit_pot_lmom(big, 0.90)
asym[["POT-Lmom u=0.90"]] <- 100 * (pot_return_levels(big, f, ex) - truth_rl) / truth_rl
asym_par <- list()
for (nm in names(COMP)) {
  th <- fit_gpd_composite(big, grid0, COMP[[nm]]$w, COMP[[nm]]$type, COMP[[nm]]$anch)
  asym_par[[nm]] <- th
  asym[[paste0(nm, " (ungrafted)")]] <-
    100 * (qgpd(1 - ex, th[1], th[2], th[3]) - truth_rl) / truth_rl
}
saveRDS(list(RL = RL, PAR = PAR, asym = asym, asym_par = asym_par,
             Tp = Tp, truth = truth_rl), "out/fits-gpd.rds")

REF <- RL[["POT-MLE u=0.90"]]          # the practitioner's standard, as reference
summ <- function(r) {
  ok <- complete.cases(r) & complete.cases(REF)
  e <- sweep(r[ok, , drop = FALSE], 2, truth_rl, "-")
  er <- sweep(REF[ok, , drop = FALSE], 2, truth_rl, "-")
  list(bias = colMeans(e), sd = apply(r[ok, , drop = FALSE], 2, sd), mse = colMeans(e^2),
       medae = apply(abs(e), 2, median), ratio = colMeans(e^2) / colMeans(er^2),
       se = apply(e^2 - er^2, 2, sd) / sqrt(nrow(e)) / colMeans(er^2),
       win = colMeans(abs(e) < abs(er)), nfail = sum(!ok))
}
S <- lapply(RL, summ)
Tshow <- c(2, 10, 20, 50, 100, 200, 500, 1000)
idx <- sapply(Tshow, function(t) which.min(abs(Tp - t)))

sink("out/gpd.txt", split = TRUE)
cat(sprintf("=== GPD study, n = %d, %d replicates ===\n", n, N_REP))
cat("Reference for all ratios: POT-MLE at the 90th percentile, hard grafted\n")
cat("onto an empirical body -- the standard practitioner's method.\n\n")
cat("--- asymptotic % bias in the return level (fit to n = 3e5) ---\n")
ab <- t(sapply(asym, function(v) v[idx])); colnames(ab) <- paste0("T=", round(Tp[idx]))
print(round(ab, 2))
cat("\n--- asymptotic composite GPD parameters (tail xi = 0.2) ---\n")
print(round(do.call(rbind, asym_par), 4))
cat("\n--- median fitted shape of the composite GPD fits ---\n")
print(round(sapply(PAR, function(m) median(m[, 3], na.rm = TRUE)), 3))
for (lab in c("ratio", "se", "medae", "win", "bias", "sd")) {
  cat(sprintf("\n--- %s ---\n", switch(lab,
    ratio = "MSE relative to POT-MLE(0.90)", se = "MC standard error of that ratio",
    medae = "median absolute error", win = "fraction closer to the truth than POT-MLE(0.90)",
    bias = "bias", sd = "sd")))
  tab <- t(sapply(S, function(s) s[[lab]][idx])); colnames(tab) <- paste0("T=", round(Tp[idx]))
  print(round(tab, 3))
}
cat("\nfailed variants:\n"); print(sapply(S, `[[`, "nfail"))
sink()

png("out/fig-gpd.png", width = 1650, height = 660, res = 133)
par(mfrow = c(1, 2), mar = c(4.4, 4.6, 3.4, 1.2))
show <- c("POT-MLE u=0.85", "POT-MLE u=0.95", "POT-Lmom u=0.90", "empirical",
          "L2, w = smoothstep(p/0.8) + graft", "L2, w = p^2 + graft",
          "L2, w = p^6 + graft", "L1, w = p^2 + graft")
cols <- c("#c0392b", "#c0392b", "#e67e22", "#7f8c8d", "#1e8449", "#52be80", "#2980b9", "#7d3c98")
ltys <- c(2, 3, 4, 3, 1, 1, 1, 1)
for (k in 1:2) {
  vals <- lapply(show, function(nm) if (k == 1) S[[nm]]$ratio else S[[nm]]$medae / S[["POT-MLE u=0.90"]]$medae)
  plot(Tp, rep(1, length(Tp)), type = "n", log = if (k == 1) "xy" else "x",
       ylim = if (k == 1) c(0.3, 20) else c(0.4, 2),
       xlab = "return period T (years)",
       ylab = if (k == 1) "MSE relative to POT-MLE(0.90)" else "median |error| relative to POT-MLE(0.90)",
       main = if (k == 1) sprintf("GPD methods (n = %d)", n) else "Median absolute error")
  abline(h = 1, col = "grey50", lwd = 1.5)
  for (i in seq_along(show)) lines(Tp, vals[[i]], col = cols[i], lwd = 2.5, lty = ltys[i])
  if (k == 1) legend("topleft", show, col = cols, lwd = 2.5, lty = ltys, bty = "n", cex = 0.6)
}
dev.off()
cat("\nwrote out/fig-gpd.png and out/gpd.txt\n")

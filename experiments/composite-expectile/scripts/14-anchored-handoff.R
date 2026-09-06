# ---------------------------------------------------------------------------
# Two follow-ups.
#
# (1) Mean anchoring. The expectile identification equation splits into a
#     tail-local partial moment and a global mean. Taking the mean from the
#     data instead of the model removes the residual contamination that no
#     weight function could: at p = 0.98 the truth's expectile sits 9.73% above
#     its GEV component's ordinarily, and 0.0001% above it once anchored.
#
# (2) Handing off later. The fitting weight and the graft's handover weight do
#     different jobs, so they are swept separately here.
#
# Output: out/fits-anchored.rds, out/anchored.txt, out/fig-anchored.png
# ---------------------------------------------------------------------------
.libPaths(c("/home/user/Rlib-graft", .libPaths()))
source("R/setup.R"); source("R/config.R"); source("R/smoothgraft.R")
suppressMessages({library(distionary); library(distplyr)})
library(parallel)

n <- N_OBS; NC <- detectCores()
Tp <- RETURN_PERIODS; ex <- 1 / Tp; truth_rl <- q_true(1 - ex)
set.seed(4242 + n)
datasets <- lapply(seq_len(N_REP), function(i) r_true(n))

# --- asymptotic targets ------------------------------------------------------
set.seed(20260906); big <- r_true(2e6)
asym <- list()
for (p0 in c(0.50, 0.90, 0.95)) {
  g <- make_level_grid(p0, N_PANEL, N_GL); wf <- make_weight(p0, min(p0 + 0.03, 0.999))
  asym[[sprintf("cee|%.2f", p0)]]      <- fit_cee(big, g, wf)
  asym[[sprintf("cee_anch|%.2f", p0)]] <- fit_cee_anchored(big, g, wf)
}
asym <- do.call(rbind, asym); colnames(asym) <- c("mu", "sigma", "xi")

# --- fits --------------------------------------------------------------------
base <- readRDS(sprintf("out/fits-n%d.rds", n))$fits
tf   <- readRDS("out/fits-tailfocus.rds")$fits
PAR <- list(mle = base$mle, lmom = base$lmom, `cee|0.50` = tf[["expectile|0.50"]])
for (p0 in c(0.50, 0.90, 0.95)) {
  g <- make_level_grid(p0, N_PANEL, N_GL); wf <- make_weight(p0, min(p0 + 0.03, 0.999))
  cat(sprintf("anchored fit p0 = %.2f ... ", p0)); t0 <- Sys.time()
  PAR[[sprintf("cee_anch|%.2f", p0)]] <- do.call(rbind, mclapply(datasets,
    function(y) fit_cee_anchored(y, g, wf, start = fit_lmom(y)), mc.cores = NC))
  cat(sprintf("%.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}

# --- grafts: sweep the handover, decoupled from the fitting weight -----------
HANDOFF <- list(c(0.90, 0.98), c(0.95, 0.99), c(0.97, 0.995), c(0.99, 0.999))
GRAFTS <- list()
for (src in c("cee|0.50", "cee_anch|0.50", "cee_anch|0.90")) {
  for (h in HANDOFF)
    GRAFTS[[sprintf("graft %s @ %.2f-%.3f", src, h[1], h[2])]] <-
      list(par = src, p0 = h[1], p1 = h[2])
}
GRAFTS[["graft mle @ 0.95-0.990"]] <- list(par = "mle", p0 = 0.95, p1 = 0.99)  # control

RL <- list()
for (nm in names(PAR))
  RL[[nm]] <- t(apply(PAR[[nm]], 1, function(th)
    if (any(is.na(th))) rep(NA_real_, length(Tp)) else qgev(1 - ex, th[1], th[2], th[3])))
RL[["empirical"]] <- do.call(rbind, mclapply(datasets, function(y)
  as.numeric(eval_quantile(dst_empirical(y), at = 1 - ex)), mc.cores = NC))
for (nm in names(GRAFTS)) {
  g <- GRAFTS[[nm]]; P <- PAR[[g$par]]
  cat(sprintf("%-34s ... ", nm)); t0 <- Sys.time()
  RL[[nm]] <- do.call(rbind, mclapply(seq_along(datasets), function(i)
    as.numeric(graft_return_levels(datasets[[i]], P[i, ], g$p0, g$p1, ex)), mc.cores = NC))
  cat(sprintf("%.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
saveRDS(list(RL = RL, PAR = PAR, asym = asym, Tp = Tp, truth = truth_rl), "out/fits-anchored.rds")

M <- RL$mle
summ <- function(r) {
  ok <- complete.cases(r) & complete.cases(M)
  e <- sweep(r[ok, , drop = FALSE], 2, truth_rl, "-")
  em <- sweep(M[ok, , drop = FALSE], 2, truth_rl, "-")
  list(bias = colMeans(e), sd = apply(r[ok, , drop = FALSE], 2, sd), mse = colMeans(e^2),
       medae = apply(abs(e), 2, median), ratio = colMeans(e^2) / colMeans(em^2),
       se = apply(e^2 - em^2, 2, sd) / sqrt(nrow(e)) / colMeans(em^2),
       win = colMeans(abs(e) < abs(em)), nfail = sum(!ok))
}
S <- lapply(RL, summ)
Tshow <- c(2, 10, 20, 50, 100, 200, 500, 1000)
idx <- sapply(Tshow, function(t) which.min(abs(Tp - t)))

sink("out/anchored.txt", split = TRUE)
cat(sprintf("=== Mean anchoring and a later handoff, n = %d, %d replicates ===\n\n", n, N_REP))
cat("--- asymptotic targets (fit to n = 2e6); true tail (0, 1, 0.2) ---\n")
print(round(asym, 4))
cat("\nasymptotic % bias in the return level\n")
ab <- t(apply(asym, 1, function(th)
  100 * (qgev(1 - ex[idx], th[1], th[2], th[3]) - truth_rl[idx]) / truth_rl[idx]))
colnames(ab) <- paste0("T=", round(Tp[idx])); print(round(ab, 2))
cat("\n--- median fitted shape (true xi = 0.2) ---\n")
print(round(sapply(PAR, function(m) median(m[, 3], na.rm = TRUE)), 3))
for (lab in c("ratio", "se", "medae", "win", "bias", "sd")) {
  cat(sprintf("\n--- %s ---\n", switch(lab, ratio = "MSE relative to GEV MLE",
    se = "MC standard error of that ratio", medae = "median absolute error",
    win = "fraction closer to the truth than the MLE", bias = "bias", sd = "sd")))
  tab <- t(sapply(S, function(s) s[[lab]][idx])); colnames(tab) <- paste0("T=", round(Tp[idx]))
  print(round(tab, 3))
}
cat("\nfailed variants:\n"); print(sapply(S, `[[`, "nfail"))
sink()

png("out/fig-anchored.png", width = 1650, height = 640, res = 133)
par(mfrow = c(1, 3), mar = c(4.4, 4.6, 3.4, 1.2))
plot(Tp, rep(1, length(Tp)), type = "n", log = "xy", ylim = c(0.4, 30),
     xlab = "return period T", ylab = "MSE relative to GEV MLE",
     main = "Mean anchoring, ungrafted")
abline(h = 1, col = "grey50", lwd = 1.5)
for (i in seq_along(c("cee|0.50", "cee_anch|0.50", "cee_anch|0.90", "cee_anch|0.95")))
  lines(Tp, S[[c("cee|0.50","cee_anch|0.50","cee_anch|0.90","cee_anch|0.95")[i]]]$ratio,
        col = c("#2980b9", "#1e8449", "#7d3c98", "#b9770e")[i], lwd = 2.6)
legend("topright", c("cee p0=0.50 (ordinary)", "anchored p0=0.50", "anchored p0=0.90",
                     "anchored p0=0.95"),
       col = c("#2980b9", "#1e8449", "#7d3c98", "#b9770e"), lwd = 2.6, bty = "n", cex = 0.7)
cols <- colorRampPalette(c("#a9cce3", "#154360"))(length(HANDOFF))
for (src in c("cee|0.50", "cee_anch|0.50")) {
  plot(Tp, rep(1, length(Tp)), type = "n", log = "xy", ylim = c(0.4, 30),
       xlab = "return period T", ylab = "MSE relative to GEV MLE",
       main = sprintf("Grafted %s:\nhandover swept", src))
  abline(h = 1, col = "grey50", lwd = 1.5)
  for (i in seq_along(HANDOFF))
    lines(Tp, S[[sprintf("graft %s @ %.2f-%.3f", src, HANDOFF[[i]][1], HANDOFF[[i]][2])]]$ratio,
          col = cols[i], lwd = 2.6)
  lines(Tp, S[["graft mle @ 0.95-0.990"]]$ratio, col = "#c0392b", lwd = 2.4, lty = 2)
  legend("topright", c(sapply(HANDOFF, function(h) sprintf("handover %.2f-%.3f", h[1], h[2])),
                       "MLE grafted (control)"),
         col = c(cols, "#c0392b"), lwd = 2.4, lty = c(rep(1, length(HANDOFF)), 2),
         bty = "n", cex = 0.68)
}
dev.off()
cat("\nwrote out/fig-anchored.png and out/anchored.txt\n")

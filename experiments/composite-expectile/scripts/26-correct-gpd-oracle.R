# ---------------------------------------------------------------------------
# 26. The correctly-specified GPD study against a real efficiency bound, and
#     the decomposition of section 13's headline gain.
#
# scripts/25 used the full-sample three-parameter GPD MLE, whose threshold
# estimate mu_hat = min(y) is a boundary parameter and therefore non-regular.
# This adds the ORACLE benchmark -- threshold fixed at the true mu, so sigma and
# xi come from an ordinary regular two-parameter MLE on all n points. That is
# the efficiency bound for the family. The two agree to 0.1% at n = 100, which
# says the non-regularity costs nothing here and the full MLE was already a
# sound benchmark.
#
# The point of all this is a confound in section 13. Its reference,
# POT-MLE(0.90), keeps ten exceedances out of a hundred observations. Measured
# against the bound it is 2.58x worse at T = 1000 even when it is perfectly
# specified. So part of section 13's result is the composite criterion and part
# is simply using more data, and the two need separating.
#
# Output: out/correct-gpd-oracle.txt, out/fig-correct-gpd.png
# ---------------------------------------------------------------------------
source("R/setup.R"); source("R/config.R")
source("R/gpd.R"); source("R/gpd_estimators.R")
library(parallel)

n <- N_OBS; NC <- detectCores()
MU <- 1; SIG <- 0.5; XI <- 0.2
Tp <- RETURN_PERIODS; ex <- 1 / Tp; truth_rl <- qgpd(1 - ex, MU, SIG, XI)
set.seed(20260907)
datasets <- lapply(seq_len(N_REP), function(i) qgpd(runif(n), MU, SIG, XI))

fit_oracle <- function(y) {
  z <- y - MU; z <- z[z > 0]
  o <- try(optim(c(log(mean(z)), 0.1), gpd_exceed_nll, z = z, method = "Nelder-Mead",
                 control = list(maxit = 2000, reltol = 1e-12)), silent = TRUE)
  if (inherits(o, "try-error") || o$value >= BIG) return(rep(NA_real_, length(ex)))
  s <- exp(o$par[1]); xi <- o$par[2]
  MU + s * (ex^(-xi) - 1) / xi
}
cat("oracle MLE ... ")
RLo <- do.call(rbind, mclapply(datasets, fit_oracle, mc.cores = NC)); cat("done\n")
prev <- readRDS("out/fits-correct-gpd-refs.rds")
RL <- c(list("ORACLE MLE (mu known)" = RLo), prev$RL)
saveRDS(list(RL = RL, Tp = Tp, truth = truth_rl), "out/fits-correct-gpd-oracle.rds")

ok0 <- complete.cases(RLo)
rel <- function(r) {
  ok <- complete.cases(r) & ok0
  e  <- sweep(r[ok, , drop = FALSE], 2, truth_rl, "-")
  eo <- sweep(RLo[ok, , drop = FALSE], 2, truth_rl, "-")
  list(ratio = colMeans(e^2) / colMeans(eo^2), bias = colMeans(e),
       win = colMeans(abs(e) < abs(eo)))
}
S <- lapply(RL, rel)
idx <- sapply(c(2, 10, 48, 107, 203, 529, 1000), function(t) which.min(abs(Tp - t)))
show <- c("ORACLE MLE (mu known)", "GPD full MLE", "POT-MLE u=0.10", "POT-MLE u=0.50",
          "POT-MLE u=0.90", "POT-Lmom u=0.90", "composite L1 + graft",
          "composite L2 + graft", "elastile a=0.5 + graft", "inv Huber k=4 + graft",
          "composite L2", "inv Huber k=4")

sink("out/correct-gpd-oracle.txt", split = TRUE)
cat(sprintf("=== Correctly specified GPD against the efficiency bound ===\n"))
cat(sprintf("iid GPD(%g, %g, %g), n = %d, %d replicates.\n\n", MU, SIG, XI, n, N_REP))
for (lab in c("ratio", "win", "bias")) {
  cat(sprintf("\n--- %s ---\n", switch(lab,
    ratio = "MSE relative to the ORACLE two-parameter MLE",
    win = "fraction closer to the truth than the oracle MLE", bias = "bias")))
  tab <- t(sapply(S[show], function(s) s[[lab]][idx]))
  colnames(tab) <- paste0("T=", round(Tp[idx])); print(round(tab, 3))
}
cat("\n--- decomposing section 13's headline gain ---\n")
cat("Section 13 reports composite L2 + graft at MSE ratio 0.209 against\n")
cat("POT-MLE(0.90) on the MISSPECIFIED truth, at T = 1000. Here, with the same\n")
cat("estimator and reference but NO misspecification at all, that ratio is:\n")
i1 <- which.min(abs(Tp - 1000))
r_correct <- { ok <- complete.cases(RL[["composite L2 + graft"]]) &
                     complete.cases(RL[["POT-MLE u=0.90"]])
  e <- sweep(RL[["composite L2 + graft"]][ok, , drop = FALSE], 2, truth_rl, "-")
  er <- sweep(RL[["POT-MLE u=0.90"]][ok, , drop = FALSE], 2, truth_rl, "-")
  (colMeans(e^2) / colMeans(er^2))[i1] }
cat(sprintf("\n   %.3f\n\n", r_correct))
lt <- log(0.209); lc <- log(r_correct)
cat(sprintf("On a log scale the total gain is %.3f, of which %.3f (%.0f%%) is present\n",
            lt, lc, 100 * lc / lt))
cat(sprintf("with a perfectly specified model and only %.3f (%.0f%%) is attributable\n",
            lt - lc, 100 * (lt - lc) / lt))
cat("to the misspecification the method exists to exploit.\n")
sink()

png("out/fig-correct-gpd.png", width = 1650, height = 620, res = 133)
par(mfrow = c(1, 2), mar = c(4.6, 4.6, 3.4, 1.2))
cols <- c("#1f5f8b", "#197a45", "#7d3c98", "#a53a2b")
nmk <- c("composite L1", "composite L2", "elastile a=0.5", "inv Huber k=4")
for (g in c(FALSE, TRUE)) {
  plot(Tp, rep(1, length(Tp)), type = "n", log = "x", ylim = c(0.6, 5),
       xlab = "return period T", ylab = "MSE relative to the oracle MLE",
       main = if (g) "grafted onto an empirical body" else "ungrafted")
  abline(h = 1, col = "grey40", lwd = 2)
  lines(Tp, S[["POT-MLE u=0.90"]]$ratio, col = "black", lwd = 3, lty = 2)
  lines(Tp, S[["POT-Lmom u=0.90"]]$ratio, col = "#e67e22", lwd = 2.4, lty = 3)
  for (i in seq_along(nmk))
    lines(Tp, S[[if (g) paste0(nmk[i], " + graft") else nmk[i]]]$ratio,
          col = cols[i], lwd = 2.4)
  legend("topright", c(nmk, "POT-MLE(0.90)", "POT-Lmom(0.90)"),
         col = c(cols, "black", "#e67e22"), lwd = 2.4,
         lty = c(rep(1, 4), 2, 3), bty = "n", cex = 0.62)
}
dev.off()
cat("\nwrote out/fig-correct-gpd.png and out/correct-gpd-oracle.txt\n")

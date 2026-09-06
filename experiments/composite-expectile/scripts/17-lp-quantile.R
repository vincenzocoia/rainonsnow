# ---------------------------------------------------------------------------
# The Lp-quantile at a = 1.5, inside the GPD study.
#
# Lp-quantiles (Chen 1996; Breckling & Chambers' M-quantiles specialised to a
# power loss; studied in the extremes literature by Daouia, Girard & Stupfler)
# interpolate the quantile (a = 1) and the expectile (a = 2) by changing the
# EXPONENT of the loss rather than mixing two losses additively. They are a
# different path between the same endpoints -- at p = 0.9 for GPD(0,1,0.2) the
# 1.5-quantile is 2.747, which is not even between the quantile (2.925) and the
# expectile (2.875).
#
# a = 1 and a = 2 are already in scripts/16 as the alpha = 0 and alpha = 1
# elastile fits, so only a = 1.5 is run here. Reference is POT-MLE(0.90) hard
# grafted, as in scripts/15 and 16.
#
# Output: out/fits-lp.rds, out/lp-quantile.txt
# ---------------------------------------------------------------------------
.libPaths(c("/home/user/Rlib-graft", .libPaths()))
source("R/setup.R"); source("R/config.R")
source("R/gpd.R"); source("R/gpd_estimators.R"); source("R/lpquantile.R")
source("R/smoothgraft.R"); source("R/graft_fast.R")
suppressMessages({library(distionary); library(distplyr)})
library(parallel)

n <- N_OBS; NC <- detectCores(); A <- 1.5
Tp <- RETURN_PERIODS; ex <- 1 / Tp; truth_rl <- q_true(1 - ex)
W <- function(p) p^6; WD <- function(p) 6 * p^5
grid0 <- make_level_grid(0, N_PANEL, N_GL)
set.seed(4242 + n)
datasets <- lapply(seq_len(N_REP), function(i) r_true(n))

cat("asymptotic target (n = 5000) ... ")
set.seed(202); big <- r_true(5000)
asym <- fit_gpd_lp_composite(big, grid0, W, A)
cat(sprintf("(%.3f, %.3f, %.3f)\n", asym[1], asym[2], asym[3]))

cat(sprintf("fitting a = %.2f on %d replicates ... ", A, N_REP)); t0 <- Sys.time()
PAR <- do.call(rbind, mclapply(datasets,
  function(y) fit_gpd_lp_composite(y, grid0, W, A), mc.cores = NC))
cat(sprintf("%.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
cat("grafting ... ")
RLg <- do.call(rbind, mclapply(seq_along(datasets), function(i) {
  th <- PAR[i, ]
  if (any(is.na(th))) return(rep(NA_real_, length(Tp)))
  as.numeric(graft_fast_return_levels(datasets[[i]], th, W, ex, "gpd", w_deriv = WD))
}, mc.cores = NC))
RLu <- t(apply(PAR, 1, function(th) if (any(is.na(th))) rep(NA_real_, length(Tp))
               else qgpd(1 - ex, th[1], th[2], th[3])))
cat("done\n")
saveRDS(list(PAR = PAR, RL_grafted = RLg, RL_ungrafted = RLu, a = A, asym = asym,
             Tp = Tp, truth = truth_rl), "out/fits-lp.rds")

el <- readRDS("out/fits-gpd-elastile.rds")
RL <- list(`POT-MLE u=0.90` = el$RL[["POT-MLE u=0.90"]],
           `POT-Lmom u=0.90` = el$RL[["POT-Lmom u=0.90"]],
           `a = 1 (quantile) + graft` = el$RL[["alpha=0.00 + graft"]],
           `a = 1.5 + graft` = RLg,
           `a = 2 (expectile) + graft` = el$RL[["alpha=1.00 + graft"]],
           `elastile alpha=0.5 + graft` = el$RL[["alpha=0.50 + graft"]])
REF <- RL[["POT-MLE u=0.90"]]
summ <- function(r) {
  ok <- complete.cases(r) & complete.cases(REF)
  e <- sweep(r[ok, , drop = FALSE], 2, truth_rl, "-")
  er <- sweep(REF[ok, , drop = FALSE], 2, truth_rl, "-")
  list(ratio = colMeans(e^2) / colMeans(er^2),
       se = apply(e^2 - er^2, 2, sd) / sqrt(nrow(e)) / colMeans(er^2),
       medae = apply(abs(e), 2, median), win = colMeans(abs(e) < abs(er)),
       nfail = sum(!ok))
}
S <- lapply(RL, summ)
Tshow <- c(20, 50, 100, 200, 500, 1000)
idx <- sapply(Tshow, function(t) which.min(abs(Tp - t)))

sink("out/lp-quantile.txt", split = TRUE)
cat(sprintf("=== Lp-quantile at a = %.1f, GPD study, n = %d, %d replicates ===\n", A, n, N_REP))
cat("w(p) = p^6 for fitting and handover. Ratios against POT-MLE(0.90) hard grafted.\n")
cat("a = 1 and a = 2 are the pure quantile and expectile composite fits.\n\n")
cat(sprintf("asymptotic fit (n = 5000): mu %.3f, sigma %.3f, xi %.3f\n", asym[1], asym[2], asym[3]))
cat(sprintf("median fitted xi at n = %d: %.3f   (fraction at the xi > 0 bound: %.3f)\n",
            n, median(PAR[, 3], na.rm = TRUE), mean(PAR[, 3] <= 1e-3, na.rm = TRUE)))
cat(sprintf("failed fits: %d of %d\n", sum(!complete.cases(PAR)), N_REP))
for (lab in c("ratio", "se", "medae", "win")) {
  cat(sprintf("\n--- %s ---\n", switch(lab, ratio = "MSE relative to POT-MLE(0.90)",
    se = "MC standard error of that ratio", medae = "median absolute error",
    win = "fraction closer to the truth than POT-MLE(0.90)")))
  tab <- t(sapply(S, function(s) s[[lab]][idx])); colnames(tab) <- paste0("T=", round(Tp[idx]))
  print(round(tab, 3))
}
cat("\n--- ungrafted, for reference ---\n")
su <- summ(RLu)
print(round(setNames(su$ratio[idx], paste0("T=", round(Tp[idx]))), 3))
sink()
cat("\nwrote out/lp-quantile.txt\n")

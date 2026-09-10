# ---------------------------------------------------------------------------
# 31. MWLE in the GPD setting, contaminated and correct.
#
# Contaminated: a three-parameter GPD fitted to the whole sample with w = p^6,
# smooth grafted, against POT-MLE(0.90) -- the protocol of section 13.
# Correct: iid GPD(0, 0.5, 0.2) with the threshold known, two parameters for
# every method, against the GPD MLE -- the protocol of section 16.
#
# Output: out/mwle-gpd.txt
# ---------------------------------------------------------------------------
.libPaths(c("/home/user/Rlib-graft", .libPaths()))
source("R/setup.R"); source("R/config.R")
source("R/gpd.R"); source("R/gpd_estimators.R"); source("R/gpd2.R")
source("R/invhuber.R"); source("R/mwle.R")
source("R/smoothgraft.R"); source("R/graft_fast.R")
suppressMessages({library(distionary); library(distplyr)})
library(parallel)

n <- N_OBS; NC <- detectCores()
Tp <- RETURN_PERIODS; ex <- 1 / Tp
W <- function(p) p^6; WD <- function(p) 6 * p^5
gr <- make_level_grid(0, N_PANEL, N_GL)
idx <- sapply(c(2, 10, 48, 203, 529, 1000), function(t) which.min(abs(Tp - t)))

report <- function(tag, RL, REF, truth_rl) {
  summ <- function(r) { ok <- complete.cases(r) & complete.cases(REF)
    e <- sweep(r[ok, , drop = FALSE], 2, truth_rl, "-")
    er <- sweep(REF[ok, , drop = FALSE], 2, truth_rl, "-")
    list(ratio = (colMeans(e^2) / colMeans(er^2))[idx],
         bias = colMeans(e)[idx], nfail = sum(!ok)) }
  S <- lapply(RL, summ)
  cat(sprintf("\n--- %s ---\n", tag))
  cat(sprintf("%-24s", "")); cat(sprintf("%8s", paste0("T=", round(Tp[idx])))); cat("\n")
  for (nm in names(RL)) { cat(sprintf("%-24s", nm)); cat(sprintf("%8.3f", S[[nm]]$ratio)); cat("\n") }
  cat("\nbias:\n")
  for (nm in names(RL)) { cat(sprintf("%-24s", nm)); cat(sprintf("%8.3f", S[[nm]]$bias)); cat("\n") }
  cat("\nfailed:\n"); print(sapply(S, `[[`, "nfail"))
}

sink("out/mwle-gpd.txt", split = TRUE)
cat(sprintf("=== MWLE in the GPD setting, n = %d, %d replicates ===\n", n, N_REP))

## ---- contaminated -------------------------------------------------------
set.seed(4242 + n)
dat <- lapply(seq_len(N_REP), function(i) r_true(n))
truth_rl <- q_true(1 - ex)
RL <- list()
RL[["POT-MLE u=0.90"]] <- do.call(rbind, mclapply(dat,
  function(y) pot_return_levels(y, fit_pot_mle(y, 0.90), ex), mc.cores = NC))
FIT3 <- list("MWLE (Fung) + graft" = function(y) fit_mwle_gpd3(y, W),
             "composite L2 + graft" = function(y) fit_gpd_composite(y, gr, W, "expectile"),
             "inv Huber k=4 + graft" = function(y) fit_gpd_onesided(y, gr, W, 4))
for (nm in names(FIT3)) {
  cat(sprintf("  %-22s ", nm)); t0 <- Sys.time()
  P <- do.call(rbind, mclapply(dat, FIT3[[nm]], mc.cores = NC))
  RL[[nm]] <- do.call(rbind, mclapply(seq_along(dat), function(i)
    as.numeric(graft_fast_return_levels(dat[[i]], P[i, ], W, ex, "gpd", w_deriv = WD)),
    mc.cores = NC))
  cat(sprintf("%.1f min   median xi %.3f\n",
      as.numeric(difftime(Sys.time(), t0, units = "mins")), median(P[, 3], na.rm = TRUE)))
}
report("CONTAMINATED: MSE relative to POT-MLE(0.90)", RL, RL[["POT-MLE u=0.90"]], truth_rl)

## ---- correctly specified ------------------------------------------------
set.seed(20260908)
dat2 <- lapply(seq_len(N_REP), function(i) qgpd(runif(n), 0, 0.5, 0.2))
truth2 <- qgpd(1 - ex, 0, 0.5, 0.2)
rl2 <- function(p) if (anyNA(p)) rep(NA_real_, length(ex)) else qgpd(1 - ex, 0, p[1], p[2])
RL2 <- list()
RL2[["GPD MLE"]]       <- do.call(rbind, mclapply(dat2, function(y) rl2(gpd2_mle(y)), mc.cores = NC))
RL2[["GPD L-moments"]] <- do.call(rbind, mclapply(dat2, function(y) rl2(gpd2_lmom(y)), mc.cores = NC))
FIT2 <- list("MWLE (Fung)"   = function(y) fit_mwle_gpd2(y, W),
             "composite L2"  = function(y) fit_gpd2_composite(y, gr, W, "expectile"),
             "inv Huber k=4" = function(y) fit_gpd2_onesided(y, gr, W, 4))
for (nm in names(FIT2)) {
  cat(sprintf("  %-22s ", nm)); t0 <- Sys.time()
  P <- do.call(rbind, mclapply(dat2, FIT2[[nm]], mc.cores = NC))
  RL2[[nm]] <- do.call(rbind, lapply(seq_len(N_REP), function(i) rl2(P[i, ])))
  RL2[[paste0(nm, " + graft")]] <- do.call(rbind, mclapply(seq_along(dat2), function(i)
    as.numeric(graft_fast_return_levels(dat2[[i]], c(0, P[i, 1], P[i, 2]), W, ex, "gpd",
                                        w_deriv = WD)), mc.cores = NC))
  cat(sprintf("%.1f min   median xi %.3f\n",
      as.numeric(difftime(Sys.time(), t0, units = "mins")), median(P[, 2], na.rm = TRUE)))
}
report("CORRECT GPD: MSE relative to the GPD MLE", RL2, RL2[["GPD MLE"]], truth2)
sink()
cat("\nwrote out/mwle-gpd.txt\n")

# ---------------------------------------------------------------------------
# 30. Fung's maximum weighted likelihood as a baseline, in the GEV setting.
#
# MWLE is an existing tail-focused estimator, so it belongs in the comparison:
# it tests whether the composite construction earns anything over simply
# tilting the likelihood towards the tail. Both truths are run -- the
# contaminated max(GEV, Normal) and an exactly correct GEV -- with the same
# weight the composite estimators use, floored away from zero as Fung requires.
#
# Output: out/mwle-gev.txt
# ---------------------------------------------------------------------------
source("R/setup.R"); source("R/config.R"); source("R/invhuber.R"); source("R/mwle.R")
library(parallel)

n <- N_OBS; NC <- detectCores()
Tp <- RETURN_PERIODS; ex <- 1 / Tp
W  <- make_weight(0.50, 0.90)                 # the weight section 5 recommends
gr <- make_level_grid(0.50, N_PANEL, N_GL)
rl <- function(p) if (anyNA(p)) rep(NA_real_, length(ex)) else qgev(1 - ex, p[1], p[2], p[3])

FITS <- list(
  "GEV MLE"        = function(y) fit_mle(y),
  "GEV L-moments"  = function(y) fit_lmom(y),
  "MWLE (Fung)"    = function(y) fit_mwle_gev(y, W),
  "composite L1"   = function(y) fit_cqe(y, gr, W),
  "composite L2"   = function(y) fit_cee(y, gr, W),
  "elastile a=0.5" = function(y) fit_elastile(y, gr, W, 0.5),
  "inv Huber k=4"  = function(y) fit_gev_onesided(y, gr, W, 4)
)

run <- function(tag, gen, truth_rl, seed) {
  set.seed(seed)
  dat <- lapply(seq_len(N_REP), function(i) gen(n))
  PAR <- list(); RL <- list()
  for (nm in names(FITS)) {
    cat(sprintf("  %-16s ", nm)); t0 <- Sys.time()
    PAR[[nm]] <- do.call(rbind, mclapply(dat, FITS[[nm]], mc.cores = NC))
    RL[[nm]] <- do.call(rbind, lapply(seq_len(N_REP), function(i) rl(PAR[[nm]][i, ])))
    cat(sprintf("%.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  }
  REF <- RL[["GEV MLE"]]
  idx <- sapply(c(2, 10, 48, 203, 529, 1000), function(t) which.min(abs(Tp - t)))
  summ <- function(r) { ok <- complete.cases(r) & complete.cases(REF)
    e <- sweep(r[ok, , drop = FALSE], 2, truth_rl, "-")
    er <- sweep(REF[ok, , drop = FALSE], 2, truth_rl, "-")
    list(ratio = (colMeans(e^2) / colMeans(er^2))[idx], bias = colMeans(e)[idx],
         xi = NA, nfail = sum(!ok)) }
  S <- lapply(RL, summ)
  cat(sprintf("\n--- %s: MSE relative to the GEV MLE ---\n", tag))
  cat(sprintf("%-17s", "")); cat(sprintf("%8s", paste0("T=", round(Tp[idx])))); cat("\n")
  for (nm in names(RL)) { cat(sprintf("%-17s", nm)); cat(sprintf("%8.3f", S[[nm]]$ratio)); cat("\n") }
  cat(sprintf("\n--- %s: bias ---\n", tag))
  for (nm in names(RL)) { cat(sprintf("%-17s", nm)); cat(sprintf("%8.3f", S[[nm]]$bias)); cat("\n") }
  cat("\nmedian fitted (mu, sigma, xi):\n")
  print(round(t(sapply(PAR, function(m) apply(m, 2, median, na.rm = TRUE))), 3))
  cat("\nfailed fits:\n"); print(sapply(S, `[[`, "nfail"))
  list(RL = RL, PAR = PAR)
}

sink("out/mwle-gev.txt", split = TRUE)
cat(sprintf("=== MWLE in the GEV setting, n = %d, %d replicates ===\n", n, N_REP))
cat("Weight: the composite study's own, floored at eps = 1e-3 for Fung's\n")
cat("positivity requirement, and mapped from p to y through Fhat.\n\n")
cat("CONTAMINATED TRUTH  max(GEV(0,1,0.2), N(1.5,0.8))\n")
a <- run("contaminated", function(k) r_true(k), q_true(1 - ex), 4242 + n)
cat("\n\nEXACTLY CORRECT GEV(0, 1, 0.2)\n")
b <- run("correct GEV", function(k) qgev(runif(k), 0, 1, 0.2),
         qgev(1 - ex, 0, 1, 0.2), 20260907)
sink()
saveRDS(list(contaminated = a, correct = b, Tp = Tp), "out/fits-mwle-gev.rds")
cat("\nwrote out/mwle-gev.txt\n")

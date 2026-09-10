# ---------------------------------------------------------------------------
# 33. The composite extremile estimator.
#
# Extremiles are L-functionals, not M-quantiles, so the composite LOSS does not
# transfer; what transfers is minimum-distance matching on the extremile
# function. That makes this a tail-weighted probability-weighted-moment
# estimator, so its natural comparison is L-moments rather than the M-estimators.
#
# The tuning parameter is not a weight but the reach: the grid runs over
# tau in [1/2, 0.5^(1/(alpha n))], so r(tau) <= alpha n. Small alpha keeps the
# empirical extremile nearly unbiased but gives up tail reach; large alpha does
# the opposite. Both truths are run, in both families.
#
# Output: out/extremile.txt
# ---------------------------------------------------------------------------
source("R/setup.R"); source("R/config.R")
source("R/gpd.R"); source("R/gpd_estimators.R"); source("R/gpd2.R")
source("R/extremile.R"); source("R/invhuber.R")
library(parallel)

n <- N_OBS; NC <- detectCores()
Tp <- RETURN_PERIODS; ex <- 1 / Tp
idx <- sapply(c(2, 10, 48, 203, 529, 1000), function(t) which.min(abs(Tp - t)))
ALPHAS <- c(0.05, 0.10, 0.20, 0.50)

blk <- function(tag, dat, truth_rl, rl_of, FITS, refname) {
  RL <- list(); PAR <- list()
  for (nm in names(FITS)) {
    cat(sprintf("  %-24s ", nm)); t0 <- Sys.time()
    PAR[[nm]] <- do.call(rbind, mclapply(dat, FITS[[nm]], mc.cores = NC))
    RL[[nm]] <- do.call(rbind, lapply(seq_along(dat), function(i) rl_of(PAR[[nm]][i, ])))
    cat(sprintf("%.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  }
  REF <- RL[[refname]]
  summ <- function(r) { ok <- complete.cases(r) & complete.cases(REF)
    e <- sweep(r[ok, , drop = FALSE], 2, truth_rl, "-")
    er <- sweep(REF[ok, , drop = FALSE], 2, truth_rl, "-")
    list(ratio = (colMeans(e^2) / colMeans(er^2))[idx], nfail = sum(!ok)) }
  S <- lapply(RL, summ)
  cat(sprintf("\n--- %s: MSE relative to %s ---\n", tag, refname))
  cat(sprintf("%-24s", "")); cat(sprintf("%8s", paste0("T=", round(Tp[idx])))); cat("\n")
  for (nm in names(RL)) { cat(sprintf("%-24s", nm)); cat(sprintf("%8.3f", S[[nm]]$ratio)); cat("\n") }
  np <- ncol(PAR[[1]])
  cat(sprintf("\nmedian fitted shape (last parameter):\n"))
  print(round(sapply(PAR, function(m) median(m[, np], na.rm = TRUE)), 3))
  cat("failed:\n"); print(sapply(S, `[[`, "nfail"))
}

sink("out/extremile.txt", split = TRUE)
cat(sprintf("=== Composite extremile estimation, n = %d, %d replicates ===\n\n", n, N_REP))

## ---- GEV, three parameters ----------------------------------------------
Wg <- make_weight(0.50, 0.90); grg <- make_level_grid(0.50, N_PANEL, N_GL)
Fg <- list("GEV MLE" = function(y) fit_mle(y), "GEV L-moments" = function(y) fit_lmom(y),
           "composite L2" = function(y) fit_cee(y, grg, Wg))
for (a in ALPHAS) Fg[[sprintf("extremile alpha=%.2f", a)]] <-
  local({ aa <- a; function(y) fit_composite_extremile(y, extremile_grid(length(y), alpha = aa),
                                                       function(t) t^6, "gev") })
rlg <- function(p) if (anyNA(p)) rep(NA_real_, length(ex)) else qgev(1 - ex, p[1], p[2], p[3])
cat("CONTAMINATED, GEV family\n"); set.seed(4242 + n)
blk("contaminated GEV", lapply(seq_len(N_REP), function(i) r_true(n)), q_true(1 - ex), rlg, Fg, "GEV MLE")
cat("\n\nCORRECT GEV(0, 1, 0.2)\n"); set.seed(20260907)
blk("correct GEV", lapply(seq_len(N_REP), function(i) qgev(runif(n), 0, 1, 0.2)),
    qgev(1 - ex, 0, 1, 0.2), rlg, Fg, "GEV MLE")

## ---- GPD, threshold known ------------------------------------------------
W6 <- function(p) p^6; grp <- make_level_grid(0, N_PANEL, N_GL)
Fp <- list("GPD MLE" = function(y) gpd2_mle(y), "GPD L-moments" = function(y) gpd2_lmom(y),
           "composite L2" = function(y) fit_gpd2_composite(y, grp, W6, "expectile"),
           "inv Huber k=4" = function(y) fit_gpd2_onesided(y, grp, W6, 4))
for (a in ALPHAS) Fp[[sprintf("extremile alpha=%.2f", a)]] <-
  local({ aa <- a; function(y) fit_composite_extremile(y, extremile_grid(length(y), alpha = aa),
                                                       W6, "gpd2") })
rlp <- function(p) if (anyNA(p)) rep(NA_real_, length(ex)) else qgpd(1 - ex, 0, p[1], p[2])
cat("\n\nCORRECT GPD(0, 0.5, 0.2), threshold known\n"); set.seed(20260908)
blk("correct GPD", lapply(seq_len(N_REP), function(i) qgpd(runif(n), 0, 0.5, 0.2)),
    qgpd(1 - ex, 0, 0.5, 0.2), rlp, Fp, "GPD MLE")
sink()
cat("\nwrote out/extremile.txt\n")

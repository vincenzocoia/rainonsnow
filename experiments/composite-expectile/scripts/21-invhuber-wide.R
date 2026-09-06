# ---------------------------------------------------------------------------
# 21. Locate the knot optimum.
#
# scripts/19 found k = 4 beating pure L2 at every return period from T = 107 to
# 1000, paired t between -5.9 and -12.5. Since c -> Inf recovers the expectile
# exactly -- the loss rho_c(u) = u^2 for u >= -c tends to u^2 everywhere -- that
# means the optimum in k is interior and sits above 4. This extends the sweep.
#
# Output: appends to out/fits-invhuber.rds as fits-invhuber-wide.rds
# ---------------------------------------------------------------------------
.libPaths(c("/home/user/Rlib-graft", .libPaths()))
source("R/setup.R"); source("R/config.R")
source("R/gpd.R"); source("R/gpd_estimators.R"); source("R/invhuber.R")
source("R/smoothgraft.R"); source("R/graft_fast.R")
suppressMessages({library(distionary); library(distplyr)})
library(parallel)

n <- N_OBS; NC <- detectCores()
Tp <- RETURN_PERIODS; ex <- 1 / Tp
set.seed(4242 + n)                                  # same datasets as scripts/16, 19
datasets <- lapply(seq_len(N_REP), function(i) r_true(n))
W <- function(p) p^6; WD <- function(p) 6 * p^5
grid0 <- make_level_grid(0, N_PANEL, N_GL)

RL <- list(); PAR <- list()
for (k in c(8, 16)) {
  nm <- sprintf("k=%.2f", k)
  cat(sprintf("%s fit ... ", nm)); t0 <- Sys.time()
  PAR[[nm]] <- do.call(rbind, mclapply(datasets,
    function(y) fit_gpd_onesided(y, grid0, W, k), mc.cores = NC))
  cat(sprintf("%.1f min, graft ... ", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  t0 <- Sys.time()
  RL[[paste0(nm, " + graft")]] <- do.call(rbind, mclapply(seq_along(datasets),
    function(i) as.numeric(graft_fast_return_levels(datasets[[i]], PAR[[nm]][i, ],
                                                    W, ex, "gpd", w_deriv = WD)),
    mc.cores = NC))
  cat(sprintf("%.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
saveRDS(list(RL = RL, PAR = PAR, Tp = Tp, truth = q_true(1 - ex)),
        "out/fits-invhuber-wide.rds")
cat("wrote out/fits-invhuber-wide.rds\n")

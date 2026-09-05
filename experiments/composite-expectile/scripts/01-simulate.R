# ---------------------------------------------------------------------------
# The Monte Carlo experiment.
#
# For each sample size, draw N_REP datasets from the contaminated truth and fit
# all four estimators to each. Datasets are generated once, up front, from a
# fixed seed, so every estimator sees exactly the same data.
#
# Output: out/fits-n<N>.rds  (N_REP x 3 matrix of fitted (mu, sigma, xi) per estimator)
# ---------------------------------------------------------------------------
source("R/setup.R")
source("R/config.R")
library(parallel)
dir.create("out", showWarnings = FALSE)

n_cores <- max(1, detectCores() - 0)
sizes <- c(50, 100, 200)

for (n in sizes) {
  set.seed(4242 + n)
  datasets <- lapply(seq_len(N_REP), function(i) r_true(n))
  cat(sprintf("n = %d: fitting %d replicates on %d cores ... ", n, N_REP, n_cores))
  t0 <- Sys.time()
  res <- mclapply(datasets, function(y) {
    lmom <- fit_lmom(y)
    list(lmom = lmom,
         mle  = fit_mle(y, start = lmom),
         cqe  = fit_cqe(y, GRID, W_FUN, start = lmom),
         cee  = fit_cee(y, GRID, W_FUN, start = lmom))
  }, mc.cores = n_cores)
  cat(sprintf("%.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))

  fits <- lapply(c("mle", "lmom", "cqe", "cee"), function(nm)
    do.call(rbind, lapply(res, `[[`, nm)))
  names(fits) <- c("mle", "lmom", "cqe", "cee")
  for (nm in names(fits)) {
    colnames(fits[[nm]]) <- c("mu", "sigma", "xi")
    nf <- sum(!complete.cases(fits[[nm]]))
    if (nf > 0) cat(sprintf("   %s: %d of %d fits failed\n", nm, nf, N_REP))
  }
  saveRDS(list(n = n, fits = fits), sprintf("out/fits-n%d.rds", n))
}
cat("done\n")

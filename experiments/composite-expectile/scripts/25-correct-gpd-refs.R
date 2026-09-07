# ---------------------------------------------------------------------------
# 25. A fair reference for the correctly-specified GPD study.
#
# scripts/24 found the composite estimators beating POT-MLE(0.90) in the far
# tail even with ZERO misspecification. That is not a win for the loss -- it is
# the reference throwing away data. At n = 100 a 0.90 threshold leaves ten
# exceedances to fit two parameters, while the composite estimators use all
# hundred points. So the comparison in section 13 is confounded: part of that
# result is the composite criterion, and part is simply using more data.
#
# The matched reference is a three-parameter GPD fitted to the whole sample.
# Its MLE is available in closed-ish form: the density (1/s)(1 + xi(x-mu)/s)^
# (-1/xi - 1) is increasing in mu for every observation, so the likelihood rises
# monotonically in mu up to min(y) and mu_hat = min(y) exactly. Conditional on
# that, sigma and xi are the ordinary two-parameter exceedance MLE. The estimator
# is non-regular -- mu_hat is superefficient and converges at rate n rather than
# sqrt(n) -- but it is the maximum likelihood estimator for the model the
# composite estimators are fitting, so it is the honest benchmark.
#
# Output: out/correct-gpd-refs.txt
# ---------------------------------------------------------------------------
source("R/setup.R"); source("R/config.R")
source("R/gpd.R"); source("R/gpd_estimators.R")
library(parallel)

n <- N_OBS; NC <- detectCores()
MU <- 1; SIG <- 0.5; XI <- 0.2
Tp <- RETURN_PERIODS; ex <- 1 / Tp
truth_rl <- qgpd(1 - ex, MU, SIG, XI)
set.seed(20260907)                                  # same datasets as scripts/24
datasets <- lapply(seq_len(N_REP), function(i) qgpd(runif(n), MU, SIG, XI))

# full-sample three-parameter GPD MLE
fit_gpd_full_mle <- function(y) {
  mu <- min(y); z <- y[y > mu] - mu
  if (length(z) < 5) return(rep(NA_real_, 3))
  o <- try(optim(c(log(mean(z)), 0.1), gpd_exceed_nll, z = z, method = "Nelder-Mead",
                 control = list(maxit = 2000, reltol = 1e-12)), silent = TRUE)
  if (inherits(o, "try-error") || o$value >= BIG) return(rep(NA_real_, 3))
  s <- exp(o$par[1]); xi <- o$par[2]
  # the fit is to exceedances of mu; zeta = (n-1)/n of the mass sits above it
  c(mu, s, xi)
}
rl_full <- function(par, y) {
  if (anyNA(par)) return(rep(NA_real_, length(ex)))
  zeta <- mean(y > par[1])
  ifelse(ex >= zeta, NA_real_,
         par[1] + par[2] * ((ex / zeta)^(-par[3]) - 1) / par[3])
}

cat("full-sample 3-parameter GPD MLE ... ")
PARF <- do.call(rbind, mclapply(datasets, fit_gpd_full_mle, mc.cores = NC))
RLF  <- do.call(rbind, lapply(seq_len(N_REP), function(i) rl_full(PARF[i, ], datasets[[i]])))
cat("done\n")
extra <- list("GPD full MLE" = RLF)
for (v in c(0.25, 0.10)) {
  cat(sprintf("POT-MLE u=%.2f ... ", v))
  extra[[sprintf("POT-MLE u=%.2f", v)]] <- do.call(rbind, mclapply(datasets,
    function(y) pot_return_levels(y, fit_pot_mle(y, v), ex), mc.cores = NC))
  cat("done\n")
}

old <- readRDS("out/fits-correct-gpd.rds")
RL <- c(old$RL, extra)
saveRDS(list(RL = RL, PAR = old$PAR, PARF = PARF, Tp = Tp, truth = truth_rl),
        "out/fits-correct-gpd-refs.rds")

REF <- RLF                                          # the efficient benchmark
summ <- function(r) {
  ok <- complete.cases(r) & complete.cases(REF)
  e  <- sweep(r[ok, , drop = FALSE], 2, truth_rl, "-")
  er <- sweep(REF[ok, , drop = FALSE], 2, truth_rl, "-")
  list(ratio = colMeans(e^2) / colMeans(er^2),
       win = colMeans(abs(e) < abs(er)), bias = colMeans(e), nfail = sum(!ok))
}
S <- lapply(RL, summ)
idx <- sapply(c(2, 10, 48, 107, 203, 529, 1000), function(t) which.min(abs(Tp - t)))
keep <- c("GPD full MLE", "POT-MLE u=0.10", "POT-MLE u=0.25", "POT-MLE u=0.50",
          "POT-MLE u=0.90", "POT-Lmom u=0.90", "composite L1", "composite L1 + graft",
          "composite L2", "composite L2 + graft", "elastile a=0.5",
          "elastile a=0.5 + graft", "inv Huber k=4", "inv Huber k=4 + graft")

sink("out/correct-gpd-refs.txt", split = TRUE)
cat(sprintf("=== Correctly specified GPD, against the FULL-SAMPLE MLE ===\n"))
cat(sprintf("iid GPD(%g, %g, %g), n = %d, %d replicates.\n", MU, SIG, XI, n, N_REP))
cat("Ratios are against the three-parameter GPD MLE fitted to the whole sample,\n")
cat("which is the efficient estimator for the model everything here is fitting.\n")
cat("Above 1 is the price paid. The POT rows show what the section 13 reference\n")
cat("costs by discarding data even when it is correctly specified.\n\n")
for (lab in c("ratio", "win", "bias")) {
  cat(sprintf("\n--- %s ---\n", switch(lab,
    ratio = "MSE relative to the full-sample GPD MLE",
    win = "fraction closer to the truth than the full-sample MLE", bias = "bias")))
  tab <- t(sapply(S[keep], function(s) s[[lab]][idx]))
  colnames(tab) <- paste0("T=", round(Tp[idx]))
  print(round(tab, 3))
}
sink()
cat("\nwrote out/correct-gpd-refs.txt\n")

# ---------------------------------------------------------------------------
# 23. The honest test: what does this cost when the family is EXACTLY right?
#
# Every result so far is on a misspecified truth -- max(GEV, Normal) -- where a
# tail-weighted criterion has misspecification to avoid. If the data really are
# iid GEV, maximum likelihood is the efficient estimator and every composite
# estimator here is throwing information away on purpose. It must lose. The
# question is by how much: that number is the insurance premium, and without it
# none of this is deployable.
#
# DGP: iid GEV(0, 1, 0.2), n = 100, 2000 replicates. Same weight, same grid,
# same return periods as the main study. The weight starts at the median, which
# is the configuration section 5 recommends.
#
# Output: out/fits-correct-gev.rds, out/correct-gev.txt, out/fig-correct-gev.png
# ---------------------------------------------------------------------------
source("R/setup.R"); source("R/config.R"); source("R/invhuber.R")
library(parallel)

n <- N_OBS; NC <- detectCores()
MU <- 0; SIG <- 1; XI <- 0.2                      # the truth, exactly a GEV
Tp <- RETURN_PERIODS; ex <- 1 / Tp
truth_rl <- qgev(1 - ex, MU, SIG, XI)

set.seed(20260907)
datasets <- lapply(seq_len(N_REP), function(i)
  qgev(runif(n), MU, SIG, XI))

W  <- make_weight(0.50, 0.90)                     # the recommended weight
gr <- make_level_grid(0.50, N_PANEL, N_GL)

rl_of <- function(par) if (anyNA(par)) rep(NA_real_, length(ex)) else
  qgev(1 - ex, par[1], par[2], par[3])

FITS <- list(
  "GEV MLE"            = function(y) fit_mle(y),
  "GEV L-moments"      = function(y) fit_lmom(y),
  "composite L1"       = function(y) fit_cqe(y, gr, W),
  "composite L2"       = function(y) fit_cee(y, gr, W),
  "elastile a=0.5"     = function(y) fit_elastile(y, gr, W, 0.5),
  "inv Huber k=4"      = function(y) fit_gev_onesided(y, gr, W, 4)
)

PAR <- list(); RL <- list()
for (nm in names(FITS)) {
  cat(sprintf("%-16s ... ", nm)); t0 <- Sys.time()
  PAR[[nm]] <- do.call(rbind, mclapply(datasets, FITS[[nm]], mc.cores = NC))
  RL[[nm]]  <- do.call(rbind, lapply(seq_len(N_REP), function(i) rl_of(PAR[[nm]][i, ])))
  cat(sprintf("%.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
saveRDS(list(PAR = PAR, RL = RL, Tp = Tp, truth = truth_rl,
             pars = c(MU, SIG, XI)), "out/fits-correct-gev.rds")

REF <- RL[["GEV MLE"]]
sq <- function(R) sweep(R, 2, truth_rl, "-")^2
summ <- function(r) {
  ok <- complete.cases(r) & complete.cases(REF)
  e <- sweep(r[ok, , drop = FALSE], 2, truth_rl, "-")
  er <- sweep(REF[ok, , drop = FALSE], 2, truth_rl, "-")
  list(ratio = colMeans(e^2) / colMeans(er^2),
       medae = apply(abs(e), 2, median), win = colMeans(abs(e) < abs(er)),
       bias = colMeans(e), nfail = sum(!ok))
}
S <- lapply(RL, summ)
idx <- sapply(c(2, 10, 50, 100, 200, 500, 1000), function(t) which.min(abs(Tp - t)))

sink("out/correct-gev.txt", split = TRUE)
cat(sprintf("=== Correctly specified GEV: iid GEV(%g, %g, %g), n = %d, %d reps ===\n",
            MU, SIG, XI, n, N_REP))
cat("The family is EXACTLY right, so the MLE is efficient and everything else\n")
cat("should lose. Ratios are against the GEV MLE. > 1 is the price paid.\n\n")
cat("--- median fitted parameters (truth 0, 1, 0.2) ---\n")
print(round(t(sapply(PAR, function(m) apply(m, 2, median, na.rm = TRUE))), 3))
for (lab in c("ratio", "medae", "win", "bias")) {
  cat(sprintf("\n--- %s ---\n", switch(lab,
    ratio = "MSE relative to the GEV MLE", medae = "median absolute error",
    win = "fraction closer to the truth than the MLE", bias = "bias")))
  tab <- t(sapply(S, function(s) s[[lab]][idx])); colnames(tab) <- paste0("T=", round(Tp[idx]))
  print(round(tab, 3))
}
cat("\nfailed fits:\n"); print(sapply(S, `[[`, "nfail"))
sink()

png("out/fig-correct-gev.png", width = 1650, height = 620, res = 133)
par(mfrow = c(1, 2), mar = c(4.6, 4.6, 3.4, 1.2))
cols <- c("black", "#e67e22", "#1f5f8b", "#197a45", "#7d3c98", "#a53a2b")
plot(Tp, rep(1, length(Tp)), type = "n", log = "x", ylim = c(0.8, 2.2),
     xlab = "return period T", ylab = "MSE relative to the GEV MLE",
     main = sprintf("Truth is exactly GEV(0, 1, %g), n = %d", XI, n))
abline(h = 1, col = "grey40", lwd = 2)
for (i in seq_along(S)) lines(Tp, S[[i]]$ratio, col = cols[i], lwd = 2.4,
                              lty = if (i == 1) 2 else 1)
legend("topleft", names(S), col = cols, lwd = 2.4,
       lty = c(2, rep(1, length(S) - 1)), bty = "n", cex = 0.7)
plot(Tp, rep(0.5, length(Tp)), type = "n", log = "x", ylim = c(0.25, 0.6),
     xlab = "return period T", ylab = "fraction closer to truth than the MLE",
     main = "Head-to-head against the MLE")
abline(h = 0.5, col = "grey40", lwd = 2)
for (i in seq_along(S)) lines(Tp, S[[i]]$win, col = cols[i], lwd = 2.4,
                              lty = if (i == 1) 2 else 1)
dev.off()
cat("\nwrote out/fig-correct-gev.png and out/correct-gev.txt\n")

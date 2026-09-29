# ---------------------------------------------------------------------------
# Does the sandwich predict the actual sampling covariance?
#
# Two settings. CORRECT: the exceedance law really is GPD(0,1,0.2), so
# theta* is the truth and the discrepancy d_p vanishes -- A should collapse to
# the outer-product term. MISSPECIFIED: exceedances of GEV(0,1,0.2) over its
# 0.75 quantile, where theta* is the pseudo-true limit and d_p is not zero, so
# the Hessian term in A is live. The threshold is held at its population value
# in both, because estimating it adds variability the theory above does not
# model.
# ---------------------------------------------------------------------------
source("penult/08-asymp.R")
library(parallel)
NC <- min(4, detectCores())
NREP <- 3000; NS <- c(500, 2000, 8000)

settings <- list(
  correct = list(P = pop_gpd(1, 0.2), u = 0),
  misspec = list(P = pop_gev(0.2),    u = pop_gev(0.2)$Qs(0.25)))

ESTS <- list(
  "composite L1, w=p^6"  = list(kind = "quantile",  w = pm(6),
                                fit = function(z, w) fit_gpd2_composite(z, GRID, w, "quantile")),
  "composite L2, w=p^6"  = list(kind = "expectile", w = pm(6),
                                fit = function(z, w) fit_gpd2_composite(z, GRID, w, "expectile")),
  "composite L2, w=p^2"  = list(kind = "expectile", w = pm(2),
                                fit = function(z, w) fit_gpd2_composite(z, GRID, w, "expectile")),
  "inv Huber k=4, w=p^6" = list(kind = "onesided",  w = pm(6),
                                fit = function(z, w) fit_gpd2_onesided(z, GRID, w, 4)))

sink("out/asymp-check.txt", split = TRUE)
cat("=== Sandwich against Monte Carlo ===\n")
cat(sprintf("%d replicates per cell; entries are sd of sqrt(n)(theta-hat - theta*)\n", NREP))
for (sn in names(settings)) {
  st <- settings[[sn]]; P <- st$P; u <- st$u
  b <- make_base(P, u); ex <- excess_on(b); Su <- P$S(u)
  cat(sprintf("\n\n############ %s : %s, u at S(u) = %.3f ############\n", sn, P$name, Su))
  for (en in names(ESTS)) {
    E <- ESTS[[en]]
    cc <- if (E$kind == "onesided") knot_c(ex, 4) else NULL
    th <- target_composite(ex, E$kind, E$w, GRID, target_lmom(ex), cc = cc)
    SW <- sandwich_cov(ex, E$kind, E$w, GRID, th, cc = cc)
    sw <- sqrt(diag(SW$V)); rho <- SW$V[1,2] / prod(sw)
    cat(sprintf("\n%s\n  theta* = (%.5f, %.5f)   max|d_p| = %.2e\n", en, th[1], th[2], SW$dp_max))
    cat(sprintf("  sandwich : sd sigma %.4f   sd xi %.4f   corr %+.3f\n", sw[1], sw[2], rho))
    for (n in NS) {
      set.seed(20260929 + n)
      M <- do.call(rbind, mclapply(seq_len(NREP), function(i) {
        z <- ex$Qc(runif(n))
        p <- try(E$fit(z, E$w), silent = TRUE)
        if (inherits(p, "try-error") || anyNA(p)) c(NA, NA) else p
      }, mc.cores = NC))
      ok <- complete.cases(M); Mo <- M[ok, , drop = FALSE]
      sm <- sqrt(n) * apply(Mo, 2, sd); rm_ <- cor(Mo)[1,2]
      bias <- sqrt(n) * (colMeans(Mo) - th)
      cat(sprintf("  n = %5d : sd sigma %.4f (%+5.1f%%)  sd xi %.4f (%+5.1f%%)  corr %+.3f  sqrt(n)*bias (%+.2f,%+.3f)  fail %d\n",
          n, sm[1], 100*(sm[1]/sw[1]-1), sm[2], 100*(sm[2]/sw[2]-1), rm_, bias[1], bias[2], sum(!ok)))
    }
  }
}
sink()
cat("\nwrote out/asymp-check.txt\n")

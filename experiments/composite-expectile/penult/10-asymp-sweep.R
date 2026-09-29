# ---------------------------------------------------------------------------
# What the asymptotic variance says, once it is available in closed form.
#
# Three questions the finite-sample studies could not answer.
#  (a) Where does each loss's variance blow up? The criterion needs E|Y| < Inf
#      for the asymmetric square, i.e. xi < 1; the sandwich needs E[Y^2] < Inf,
#      i.e. xi < 1/2. The pinball loss needs neither. So the two should separate
#      sharply as xi approaches 1/2, and that is a prediction, not a fit.
#  (b) At matched weight, which loss is actually more efficient? The report's
#      headline is that L2 beats L1 two- to six-fold at n = 100. The asymptotic
#      variance is a different question with, it turns out, a different answer.
#  (c) What does the weight cost? Every effective threshold in section 21 is
#      bought with variance; this prices it.
# ---------------------------------------------------------------------------
source("penult/08-asymp.R")

sd_of <- function(ex, kind, w, cc = NULL) {
  th <- target_composite(ex, kind, w, GRID, target_lmom(ex), cc = cc)
  SW <- try(sandwich_cov(ex, kind, w, GRID, th, cc = cc), silent = TRUE)
  if (inherits(SW, "try-error")) return(c(NA, NA, NA, NA))
  c(th[2], sqrt(diag(SW$V)), SW$V[1,2]/prod(sqrt(diag(SW$V))))
}

sink("out/asymp-sweep.txt", split = TRUE)
cat("=== (a) The moment conditions, priced ===\n")
cat("Asymptotic sd of sqrt(n)*xi-hat on an exactly correct GPD(0,1,xi), w = p^6.\n")
cat("The asymmetric square needs xi < 1/2 for a finite sandwich; the pinball\n")
cat("loss needs nothing. Ratio is L2 / L1.\n\n")
XIS <- c(0.05, 0.1, 0.2, 0.3, 0.4, 0.45, 0.47, 0.49)
cat(sprintf("%8s %10s %10s %8s\n", "xi", "L1", "L2", "ratio"))
for (xi in XIS) {
  ex <- excess_on(make_base(pop_gpd(1, xi), 0))
  a <- sd_of(ex, "quantile", pm(6)); b <- sd_of(ex, "expectile", pm(6))
  cat(sprintf("%8.2f %10.3f %10.3f %8.2f\n", xi, a[3], b[3], b[3]/a[3]))
}

cat("\n\n=== (b) and (c) The weight, priced, on GPD(0,1,0.2) ===\n")
cat("Asymptotic sd of sqrt(n)*xi-hat, against the effective threshold each\n")
cat("weight buys (section 21, same population).\n\n")
ex <- excess_on(make_base(pop_gpd(1, 0.2), 0))
bG <- make_base(pop_gev(0.2), pop_gev(0.2)$Qs(0.25))     # for R_e, GEV xi=0.2
exG <- excess_on(bG)
cat(sprintf("%-22s %9s %9s %9s\n", "estimator", "sd(xi)", "rel MLE", "R_e (GEV)"))
ref <- NA
rows <- list()
for (m in c(0, 1, 2, 4, 6, 10, 20)) {
  for (kind in c("quantile", "expectile")) {
    s <- sd_of(ex, kind, pm(m))
    thG <- target_composite(exG, kind, pm(m), GRID, target_lmom(exG))
    Re <- effective_threshold(bG, thG[2], thG[1])$R
    rows[[length(rows)+1]] <- data.frame(
      est = sprintf("%s, w = p^%d", ifelse(kind=="quantile","L1","L2"), m),
      sd = s[3], Re = Re)
  }
}
## MLE reference: the GPD MLE's asymptotic sd for xi is sqrt((1+xi)^2)
mle_sd <- 1 + 0.2
D <- do.call(rbind, rows)
for (i in seq_len(nrow(D)))
  cat(sprintf("%-22s %9.3f %9.2f %9.1f\n", D$est[i], D$sd[i], D$sd[i]/mle_sd, D$Re[i]))
cat(sprintf("%-22s %9.3f %9.2f %9.1f\n", "GPD MLE (exact)", mle_sd, 1.00,
            effective_threshold(bG, target_mle(exG)[2], target_mle(exG)[1])$R))

cat("\n\n=== the one-sided loss, knot sweep, GPD(0,1,0.2), w = p^6 ===\n")
cat(sprintf("%8s %10s %10s\n", "k", "sd(xi)", "rel L2"))
l2 <- sd_of(ex, "expectile", pm(6))[3]
for (k in c(0.5, 1, 2, 4, 8, 16, 64)) {
  s <- sd_of(ex, "onesided", pm(6), cc = knot_c(ex, k))
  cat(sprintf("%8.1f %10.3f %10.3f\n", k, s[3], s[3]/l2))
}
sink()
cat("\nwrote out/asymp-sweep.txt\n")

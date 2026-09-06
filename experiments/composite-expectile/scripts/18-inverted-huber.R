# ---------------------------------------------------------------------------
# 18. Inverted-Huber M-quantiles: does capping the downside remove the mean
#     anchoring?
#
# Section 3 of the report established the obstacle: the expectile identification
# equation carries the model's global mean, so contamination in the body never
# fully leaves the fitted tail, and no weight function repairs it. Section 10
# showed the direct fix -- take the mean from the data instead of the model --
# works asymptotically and loses at n = 100 because the sample mean is too noisy.
#
# The inverted-Huber psi offers a structural fix instead of a plug-in one. This
# script measures, for both inversions, how much of the body contamination
# survives into the functional at each level -- the same diagnostic that gave
# 9.73% for the expectile and 0% for the quantile.
#
# Output: out/invhuber-population.txt, out/fig-invhuber-population.png
# ---------------------------------------------------------------------------
source("R/setup.R")
source("R/invhuber.R")

tr <- TRUTH
m_tr <- mean_true(tr)
phi_tr <- function(x) vapply(x, function(xx) partial_moment_true(xx, tr), numeric(1))
F_tr   <- function(x) p_true(x, tr)

# the GEV component alone: the "uncontaminated" comparison
phi_g <- function(x) gev_partial_moment(x, tr$mu, tr$sigma, tr$xi)
F_g    <- function(x) pgev(x, tr$mu, tr$sigma, tr$xi)
m_g    <- gev_mean(tr$mu, tr$sigma, tr$xi)

# surviving contamination: how far the truth's functional sits above the GEV
# component's, as a fraction. 0% means the body has been fully shed.
surv <- function(t_truth, t_gev) 100 * (t_truth / t_gev - 1)

PS <- c(0.90, 0.95, 0.98, 0.99)
CS <- c(0.1, 0.25, 0.5, 1, 2, 4, 8, 16, 32)
LO <- 0.01; HI <- 400

qs <- q_true(PS, tr); qg <- qgev(PS, tr$mu, tr$sigma, tr$xi)
es <- expectile_true(PS, tr); eg <- gev_expectile(PS, tr$mu, tr$sigma, tr$xi)

sink("out/invhuber-population.txt", split = TRUE)
cat("=== Inverted-Huber M-quantiles: surviving body contamination ===\n\n")
cat("Truth: max(GEV(0,1,0.2), N(1.5,0.8)). Entries are 100*(T_truth/T_gev - 1):\n")
cat("the percentage by which the contaminated truth's functional exceeds its own\n")
cat("GEV component's. 0 means the body has been shed entirely.\n\n")

cat("reference:\n")
cat(sprintf("  %-24s", "quantile"));  cat(sprintf("%9.4f", surv(qs, qg))); cat("\n")
cat(sprintf("  %-24s", "expectile")); cat(sprintf("%9.4f", surv(es, eg))); cat("\n")
cat(sprintf("  %-24s", "p ="));       cat(sprintf("%9.2f", PS)); cat("\n\n")

cat("(A) symmetric inverted Huber, psi(u) = sign(u) max(|u|, c):\n")
A <- matrix(NA_real_, length(CS), length(PS))
for (i in seq_along(CS)) {
  ts <- invhuber_general(PS, CS[i], phi_tr, F_tr, m_tr, LO, HI)
  tg <- invhuber_general(PS, CS[i], phi_g,  F_g,  m_g,  LO, HI)
  A[i, ] <- surv(ts, tg)
  cat(sprintf("  c = %-6.2f            ", CS[i])); cat(sprintf("%9.4f", A[i, ])); cat("\n")
}
cat("\n(B) one-sided, psi(u) = max(u, -c):\n")
B <- matrix(NA_real_, length(CS), length(PS))
for (i in seq_along(CS)) {
  ts <- onesided_general(PS, CS[i], phi_tr, LO, HI)
  tg <- onesided_general(PS, CS[i], phi_g,  LO, HI)
  B[i, ] <- surv(ts, tg)
  cat(sprintf("  c = %-6.2f            ", CS[i])); cat(sprintf("%9.4f", B[i, ])); cat("\n")
}
sink()

png("out/fig-invhuber-population.png", width = 1650, height = 640, res = 133)
par(mfrow = c(1, 2), mar = c(4.4, 4.4, 3, 1))
for (k in 1:2) {
  M <- if (k == 1) A else B
  ttl <- if (k == 1) "(A) symmetric  psi(u) = sign(u) max(|u|, c)" else
                     "(B) one-sided  psi(u) = max(u, -c)"
  matplot(CS, M, type = "b", log = "x", pch = 16, lty = 1,
          col = c("#1f5f8b", "#197a45", "#a53a2b", "#8a6d1f"),
          xlab = "knot c", ylab = "surviving contamination (%)", main = ttl)
  abline(h = 0, col = "grey60", lty = 2)
  for (j in seq_along(PS))
    abline(h = surv(es, eg)[j], col = c("#1f5f8b","#197a45","#a53a2b","#8a6d1f")[j], lty = 3)
  legend(if (k == 1) "topleft" else "topright",
         sprintf("p = %.2f", PS), col = c("#1f5f8b","#197a45","#a53a2b","#8a6d1f"),
         lty = 1, pch = 16, bty = "n", cex = 0.85)
}
dev.off()
cat("\nwrote out/fig-invhuber-population.png and out/invhuber-population.txt\n")

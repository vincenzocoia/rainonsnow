# ---------------------------------------------------------------------------
# 22. The knot sweep, drawn. Combines scripts/19 (k = 0.25 .. 4) and scripts/21
#     (k = 8, 16) with the elastile study's pure L2 and alpha = 0.5 rows, all on
#     the same 2000 datasets.
#
# Output: out/fig-invhuber-gpd.png, out/invhuber-sweep.txt
# ---------------------------------------------------------------------------
# (no setup needed: this script only reads saved fits)
ih <- readRDS("out/fits-invhuber.rds"); iw <- readRDS("out/fits-invhuber-wide.rds")
el <- readRDS("out/fits-gpd-elastile.rds")
Tp <- ih$Tp; tru <- ih$truth; REF <- ih$RL[["POT-MLE u=0.90"]]
KS <- c(0.25, 0.5, 1, 2, 4, 8, 16)
RLk <- c(ih$RL[sprintf("k=%.2f + graft", c(0.25,0.5,1,2,4))],
         iw$RL[sprintf("k=%.2f + graft", c(8,16))])
L2 <- el$RL[["alpha=1.00 + graft"]]; EL <- el$RL[["alpha=0.50 + graft"]]
sq <- function(R) sweep(R, 2, tru, "-")^2
ratio <- function(R) { ok <- complete.cases(R) & complete.cases(REF)
  colMeans(sq(R)[ok, , drop = FALSE]) / colMeans(sq(REF)[ok, , drop = FALSE]) }
Rk <- lapply(RLk, ratio); rL2 <- ratio(L2); rEL <- ratio(EL)

png("out/fig-invhuber-gpd.png", width = 1650, height = 640, res = 133)
par(mfrow = c(1, 2), mar = c(4.6, 4.6, 3.4, 1.2))
cols <- colorRampPalette(c("#a53a2b", "#7d3c98", "#197a45"))(length(KS))
plot(Tp, rep(1, length(Tp)), type = "n", log = "xy", ylim = c(0.12, 6),
     xlab = "return period T", ylab = "MSE relative to POT-MLE(0.90)",
     main = "One-sided inverted Huber + graft (n = 100)")
abline(h = 1, col = "grey50", lwd = 1.5)
for (i in seq_along(KS)) lines(Tp, Rk[[i]], col = cols[i], lwd = 2.2)
lines(Tp, rL2, col = "black", lwd = 3, lty = 2)
lines(Tp, rEL, col = "#e67e22", lwd = 3, lty = 3)
legend("bottomleft", c(sprintf("k = %g", KS), "pure L2 (k = Inf)", "elastile 0.5"),
       col = c(cols, "black", "#e67e22"), lwd = 2.2,
       lty = c(rep(1, length(KS)), 2, 3), bty = "n", cex = 0.62)

# the U in k, at two long return periods
kk <- c(which.min(abs(Tp - 529)), which.min(abs(Tp - 1000)))
plot(NA, xlim = c(0.9, 40), ylim = c(0.15, 0.46), log = "x",
     xlab = "knot k   (c = k x IQR;  k -> Inf is pure L2)",
     ylab = "MSE relative to POT-MLE(0.90)", main = "The optimum is interior")
for (j in seq_along(kk)) {
  cl <- c("#1f5f8b", "#a53a2b")[j]
  v <- sapply(Rk, function(r) r[kk[j]])
  lines(KS, v, type = "b", pch = 19, lwd = 2.6, col = cl)
  abline(h = rL2[kk[j]], col = cl, lwd = 1.8, lty = 2)
  points(4, v[KS == 4], pch = 21, cex = 2.2, lwd = 2.4, col = cl, bg = NA)
  text(30, rL2[kk[j]] + 0.012, sprintf("pure L2, T = %.0f", Tp[kk[j]]),
       cex = 0.6, col = cl, adj = 1)
}
legend("topright", sprintf("T = %.0f", Tp[kk]), col = c("#1f5f8b", "#a53a2b"),
       lwd = 2.6, pch = 19, bty = "n", cex = 0.72)
dev.off()

sink("out/invhuber-sweep.txt", split = TRUE)
cat("=== Inverted-Huber knot sweep, GPD + graft, n = 100, 2000 replicates ===\n")
cat("psi(u) = max(u, -c), c = k * IQR(y). k -> Inf recovers pure L2 exactly.\n\n")
idx <- sapply(c(19,48,107,203,529,1000), function(t) which.min(abs(Tp - t)))
tab <- rbind(do.call(rbind, lapply(Rk, function(r) r[idx])),
             "pure L2 (k=Inf)" = rL2[idx], "elastile 0.5" = rEL[idx])
rownames(tab)[seq_along(KS)] <- sprintf("k = %g", KS)
colnames(tab) <- paste0("T=", round(Tp[idx]))
cat("--- MSE relative to POT-MLE(0.90) ---\n"); print(round(tab, 3))
cat("\n--- median fitted shape (true GEV tail xi = 0.2) ---\n")
print(round(c(sapply(ih$PAR, function(m) median(m[,3], na.rm=TRUE)),
              sapply(iw$PAR, function(m) median(m[,3], na.rm=TRUE))), 3))
cat("\n--- paired MSE difference against pure L2, with t ---\n")
cat("negative = inverted Huber better\n")
for (i in seq_along(KS)) {
  A <- sq(RLk[[i]])[, idx, drop=FALSE]; B <- sq(L2)[, idx, drop=FALSE]
  ok <- complete.cases(A) & complete.cases(B); D <- A[ok,,drop=FALSE] - B[ok,,drop=FALSE]
  m <- colMeans(D); se <- apply(D, 2, sd)/sqrt(sum(ok))
  cat(sprintf("k = %-5g", KS[i])); cat(sprintf(" %7.2f[%6.1f]", m, m/se)); cat("\n")
}
sink()
cat("\nwrote out/fig-invhuber-gpd.png and out/invhuber-sweep.txt\n")

# Re-draw section 19's figure from the saved fits (script 32's plot block only).
source("R/setup.R"); source("R/config.R")
z <- readRDS("out/fits-hard-vs-smooth.rds")
RL <- z$RL; Tp <- z$Tp; truth_rl <- z$truth
REF <- RL[["POT-MLE u=0.90"]]; ok <- complete.cases(REF)
VS <- c(0.50, 0.70, 0.80, 0.90, 0.95)
full <- function(r) { o <- complete.cases(r) & ok
  e <- sweep(r[o, , drop = FALSE], 2, truth_rl, "-")
  er <- sweep(REF[o, , drop = FALSE], 2, truth_rl, "-")
  colMeans(e^2) / colMeans(er^2) }
png("out/fig-hard-vs-smooth.png", width = 1700, height = 700, res = 133)
par(mfrow = c(1, 2), mar = c(4.6, 4.8, 3.4, 1.2))
cols <- colorRampPalette(c("#e8b96a", "#a53a2b"))(length(VS))
for (nm in c("composite L2", "inv Huber k=4")) {
  plot(Tp, rep(1, length(Tp)), type = "n", log = "xy", ylim = c(0.13, 13),
       xlab = "return period T", ylab = "MSE relative to POT-MLE(0.90)", main = nm)
  abline(h = 1, col = "grey55", lwd = 1.5)
  lines(Tp, full(RL[[paste0(nm, " | ungrafted")]]), col = "#1f5f8b", lwd = 3, lty = 3)
  for (i in seq_along(VS))
    lines(Tp, full(RL[[sprintf("%s | hard v=%.2f", nm, VS[i])]]), col = cols[i], lwd = 2.2)
  lines(Tp, full(RL[[paste0(nm, " | SMOOTH")]]), col = "#197a45", lwd = 3.4)
  legend("topright", c("ungrafted", sprintf("hard, v = %.2f", VS), "smooth graft"),
         col = c("#1f5f8b", cols, "#197a45"), lwd = c(3, rep(2.2, length(VS)), 3.4),
         lty = c(3, rep(1, length(VS)), 1), bty = "n", cex = 0.68)
}
dev.off()
cat("wrote out/fig-hard-vs-smooth.png\n")

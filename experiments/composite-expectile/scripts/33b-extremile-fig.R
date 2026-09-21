# Figure for section 20: the extremile alpha-sweep, full return-period curves.
source("R/setup.R"); source("R/config.R")
source("R/gpd.R"); source("R/gpd_estimators.R"); source("R/gpd2.R"); source("R/extremile.R")
library(parallel)
n <- N_OBS; NC <- detectCores(); Tp <- RETURN_PERIODS; ex <- 1 / Tp
ALPHAS <- c(0.02, 0.05, 0.10, 0.20, 0.50)
W6 <- function(t) t^6

curves <- function(dat, truth_rl, fits, rl_of, refname) {
  RL <- lapply(fits, function(f) {
    P <- do.call(rbind, mclapply(dat, f, mc.cores = NC))
    do.call(rbind, lapply(seq_len(nrow(P)), function(i) rl_of(P[i, ]))) })
  names(RL) <- names(fits); REF <- RL[[refname]]; ok <- complete.cases(REF)
  lapply(RL, function(r) { o <- complete.cases(r) & ok
    e <- sweep(r[o, , drop = FALSE], 2, truth_rl, "-")
    er <- sweep(REF[o, , drop = FALSE], 2, truth_rl, "-")
    colMeans(e^2) / colMeans(er^2) })
}

## correct GPD(0, 0.5, 0.2), threshold known
set.seed(20260908)
d2 <- lapply(seq_len(N_REP), function(i) qgpd(runif(n), 0, 0.5, 0.2))
t2 <- qgpd(1 - ex, 0, 0.5, 0.2)
f2 <- list("GPD MLE" = gpd2_mle, "GPD L-moments" = gpd2_lmom)
for (a in ALPHAS) f2[[sprintf("a=%.2f", a)]] <- local({ aa <- a; function(y)
  fit_composite_extremile(y, extremile_grid(length(y), alpha = aa), W6, "gpd2") })
C2 <- curves(d2, t2, f2, function(p) if (anyNA(p)) rep(NA_real_, length(ex)) else
  qgpd(1 - ex, 0, p[1], p[2]), "GPD MLE")

## contaminated GEV
set.seed(4242 + n)
dg <- lapply(seq_len(N_REP), function(i) r_true(n))
tg <- q_true(1 - ex)
fg <- list("GEV MLE" = function(y) fit_mle(y), "GEV L-moments" = function(y) fit_lmom(y))
for (a in ALPHAS) fg[[sprintf("a=%.2f", a)]] <- local({ aa <- a; function(y)
  fit_composite_extremile(y, extremile_grid(length(y), alpha = aa), W6, "gev") })
CG <- curves(dg, tg, fg, function(p) if (anyNA(p)) rep(NA_real_, length(ex)) else
  qgev(1 - ex, p[1], p[2], p[3]), "GEV MLE")

png("out/fig-extremile.png", width = 1700, height = 700, res = 133)
par(mfrow = c(1, 2), mar = c(4.6, 4.8, 3.4, 1.2))
cols <- colorRampPalette(c("#cfe3f0", "#12496e"))(length(ALPHAS))
pan <- function(C, ref, main, ylim) {
  plot(Tp, rep(1, length(Tp)), type = "n", log = "xy", ylim = ylim,
       xlab = "return period T", ylab = paste("MSE relative to", ref), main = main)
  abline(h = 1, col = "grey55", lwd = 1.5)
  for (i in seq_along(ALPHAS)) lines(Tp, C[[sprintf("a=%.2f", ALPHAS[i])]], col = cols[i], lwd = 2.6)
  lines(Tp, C[[2]], col = "#d1791f", lwd = 2.6, lty = 2)
  legend("topleft", c(sprintf("extremile, alpha = %.2f", ALPHAS), "L-moments"),
         col = c(cols, "#d1791f"), lwd = 2.6, lty = c(rep(1, length(ALPHAS)), 2),
         bty = "n", cex = 0.7)
}
pan(C2, "the GPD MLE", "correct GPD(0, 0.5, 0.2)", c(0.6, 3))
pan(CG, "the GEV MLE", "contaminated truth, GEV family", c(0.4, 3))
dev.off()
cat("wrote out/fig-extremile.png\n")

idx <- sapply(c(2, 10, 48, 203, 529, 1000), function(t) which.min(abs(Tp - t)))
sink("out/extremile-fig.txt", split = TRUE)
cat("alpha sweep including 0.02, MSE ratios at selected T\n\n")
cat("--- correct GPD, relative to the GPD MLE ---\n")
cat(sprintf("%-16s", "")); cat(sprintf("%8s", paste0("T=", round(Tp[idx])))); cat("\n")
for (nm in names(C2)) { cat(sprintf("%-16s", nm)); cat(sprintf("%8.3f", C2[[nm]][idx])); cat("\n") }
cat("\n--- contaminated GEV, relative to the GEV MLE ---\n")
cat(sprintf("%-16s", "")); cat(sprintf("%8s", paste0("T=", round(Tp[idx])))); cat("\n")
for (nm in names(CG)) { cat(sprintf("%-16s", nm)); cat(sprintf("%8.3f", CG[[nm]][idx])); cat("\n") }
sink()

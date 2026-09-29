# Tables and figures for the penultimate study.
source("penult/04-setup.R")
z <- readRDS("out/penult.rds"); TAB <- z$TAB; RES <- z$RES

# r' is non-monotone for the LP3 pair with beta = 0.45 -- it rises from 0.459 at
# R = 4 to about 0.480 near R = 20 and decays back over the next eight decades --
# so r'(x) = xi* has two roots and the effective threshold is not identified for
# any estimator whose target lands inside that window. Those two populations are
# reported separately rather than quoted.
MONO <- c("GEV xi=0.00", "GEV xi=0.20", "GEV xi=0.45", "LP3 a=0.6 b=0.20", "LP3 a=1.6 b=0.20")
NONMONO <- setdiff(unique(TAB$pop), MONO)
EST <- unique(TAB$est)

wide <- function(v, pops) {
  m <- matrix(NA_real_, length(EST), length(pops), dimnames = list(EST, pops))
  for (i in seq_len(nrow(TAB))) if (TAB$pop[i] %in% pops)
    m[TAB$est[i], TAB$pop[i]] <- TAB[[v]][i]
  m
}

sink("out/penult-report.txt", split = TRUE)
cat("=================================================================\n")
cat(" EFFECTIVE THRESHOLD AND KERNEL, composite family vs references\n")
cat(" Population level (no Monte Carlo); all estimators fit above R = 4.\n")
cat("=================================================================\n\n")
cat("PIPELINE CHECK against the smooth-graft project's independent numbers\n")
mle <- TAB[TAB$est == "MLE" & TAB$pop %in% MONO, ]
cat(sprintf("  MLE effective threshold here: R = %.1f to %.1f   (their reference 13-16)\n",
            min(mle$R_e), max(mle$R_e)))
cat(sprintf("  MLE kernel 10%% point:         R = %.1f to %.1f   (their reference 6.5)\n",
            min(mle$q10), max(mle$q10)))
cat(sprintf("  MLE kernel 90%% point:         R = %.0f to %.0f   (their reference 142;\n",
            min(mle$q90), max(mle$q90)))
cat(sprintf("                                 LP3(0.6,0.20) gives %.0f)\n",
            mle$q90[mle$pop == "LP3 a=0.6 b=0.20"]))
cat(sprintf("  scale check |sigma_u/sigma* - 1|: max %.1f%%      (their tolerance ~2.5%%)\n\n",
            100 * max(abs(TAB$sig_err[TAB$pop %in% MONO]), na.rm = TRUE)))

cat("-----------------------------------------------------------------\n")
cat("1. EFFECTIVE THRESHOLD (return period at which r' equals the target shape)\n")
cat("-----------------------------------------------------------------\n")
m <- wide("R_e", MONO)
print(round(m, 1))
cat("\nreference points supplied from the smooth-graft project:\n")
cat("  MLE 13-16   WCL (linear weight) 17-22   AD2R 20-24   ADR ~11\n")

cat("\n-----------------------------------------------------------------\n")
cat("2. KERNEL SPREAD, as return periods (GEV xi = 0.20)\n")
cat("-----------------------------------------------------------------\n")
g <- TAB[TAB$pop == "GEV xi=0.20", ]
g$span <- g$q90 / g$q10
o <- order(g$R_e)
print(format(data.frame(estimator = g$est[o], xi = round(g$xi[o], 4),
                        R_e = round(g$R_e[o], 1),
                        K10 = round(g$q10[o], 1), K25 = round(g$q25[o], 1),
                        K50 = round(g$q50[o], 1), K75 = round(g$q75[o], 1),
                        K90 = round(g$q90[o], 0),
                        span = round(g$span[o], 0), intK = round(g$intK[o], 3)),
             digits = 4), row.names = FALSE)

cat("\n-----------------------------------------------------------------\n")
cat("3. THE TWO POPULATIONS WHERE THE READING IS NOT AVAILABLE\n")
cat("-----------------------------------------------------------------\n")
cat("LP3 with beta = 0.45 has non-monotone r': it rises to about 0.480 near\n")
cat("R = 20 and decays back towards 0.45 over the following decades, so\n")
cat("r'(x) = xi* has two roots for any target in that window. 51 of the 52\n")
cat("estimator-population cells there report multiple crossings. The kernel is\n")
cat("still well defined; the single number is not.\n\n")
print(round(wide("R_e", NONMONO), 1))
sink()

## ---- figure ---------------------------------------------------------------
png("out/fig-penult.png", width = 1750, height = 780, res = 132)
par(mfrow = c(1, 2), mar = c(4.4, 4.6, 3.2, 1.0))

## left: effective threshold against tail weighting
pops <- MONO; cols <- c("#1f5f8b", "#2d6a4a", "#9c4a2f", "#7a5ea8", "#b08d20")
ms <- c(0, 1, 2, 4, 6, 10, 20)
plot(range(ms + 1), c(8, 130), type = "n", log = "xy", xaxt = "n",
     xlab = expression(paste("weight exponent m in w(p) = ", p^m, "   (1 = flat)")),
     ylab = "effective threshold, return period",
     main = "Where the weight puts the estimator")
axis(1, at = ms + 1, labels = ms)
rect(1, 17, 21, 22, col = "#00000010", border = NA)
rect(1, 20, 21, 24, col = "#00000010", border = NA)
text(16, 19.5, "WCL / AD2R", cex = 0.62, col = "grey30")
for (i in seq_along(pops)) {
  v <- vapply(ms, function(mm) TAB$R_e[TAB$pop == pops[i] &
        TAB$est == sprintf("L2 (expectile) p^%d", mm)], numeric(1))
  lines(ms + 1, v, col = cols[i], lwd = 2.6)
  points(1, TAB$R_e[TAB$pop == pops[i] & TAB$est == "MLE"], col = cols[i], pch = 19, cex = 1.1)
}
legend("topleft", c(pops, "filled circle = MLE"), col = c(cols, NA), lwd = c(rep(2.6, 5), NA),
       bty = "n", cex = 0.66)

## right: the kernels themselves
kk <- RES[["GEV xi=0.20"]]$kern
show <- c("MLE", "L2 (expectile) p^0", "L2 (expectile) p^6", "inv Huber k=4 p^6",
          "MWLE p^6", "extremile rmax=50")
cl <- c("#333333", "#8fb8d4", "#1f5f8b", "#9c4a2f", "#2d6a4a", "#b08d20")
xr <- c(4, 3e4)
plot(xr, c(0, 1.05), type = "n", log = "x", xlab = "return period",
     ylab = "kernel density (per unit log R, rescaled)",
     main = "Which part of the tail each one reads")
for (i in seq_along(show)) {
  k <- kk[[show[i]]]; if (is.null(k)) next
  d <- pmax(k$dens_e, 0)                 # density per unit e = per unit log R
  lines(k$R, d / max(d), col = cl[i], lwd = 2.6)
}
legend("topright", show, col = cl, lwd = 2.6, bty = "n", cex = 0.66)
dev.off()
cat("\nwrote out/penult-report.txt and out/fig-penult.png\n")

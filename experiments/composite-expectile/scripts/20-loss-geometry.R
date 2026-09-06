# ---------------------------------------------------------------------------
# 20. The geometry of mixing losses, and why the elastic-net picture does not
#     transfer.
#
# The elastic net's diamond-and-circle diagram is a picture of a PENALTY's level
# set in coefficient space: the corners of the L1 diamond are where the solution
# lands on an axis, which is sparsity. The elastile mixes L1 and L2 on the
# RESIDUAL, not on the coefficient, so that diagram is about a different object.
#
# Left panel: the analogue does exist, in residual space. The set
# {r : sum_i rho(r_i) <= 1} is a diamond for L1, a disc for L2, and the same
# rounded diamond as the elastic net in between. What its corner encodes is not
# sparse coefficients but sparse RESIDUALS -- the fit passing exactly through
# data points, which is the classical "exact fit" property of L1 regression.
#
# Right panel: the picture that actually matters for a loss is the influence
# function psi = rho'/2. Robust statistics selects on one property of it --
# boundedness, which is what gives finite gross-error sensitivity -- and the
# additive L1/L2 mixture fails that test for every alpha > 0. Huber's splice
# passes it. That is why the field took the splice and not the mixture.
#
# Output: out/fig-loss-geometry.png
# ---------------------------------------------------------------------------
png("out/fig-loss-geometry.png", width = 1650, height = 700, res = 133)
par(mfrow = c(1, 2), mar = c(4.6, 4.6, 3.4, 1.2))

# --- residual-space unit balls ---------------------------------------------
rho_mix <- function(u, a) a * u^2 + (1 - a) * abs(u)
AL <- c(0, 0.25, 0.5, 0.75, 1)
cols <- colorRampPalette(c("#1f5f8b", "#7d3c98", "#197a45"))(length(AL))
th <- seq(0, 2 * pi, length.out = 1200)
plot(NA, xlim = c(-1.15, 1.15), ylim = c(-1.15, 1.15), asp = 1,
     xlab = expression(r[1]), ylab = expression(r[2]),
     main = "Loss balls in RESIDUAL space")
abline(h = 0, v = 0, col = "grey85")
for (i in seq_along(AL)) {
  a <- AL[i]
  # radius along each direction: solve rho_mix(t*cos) + rho_mix(t*sin) = 1
  rad <- vapply(th, function(t0) {
    d <- c(cos(t0), sin(t0))
    f <- function(r) sum(rho_mix(r * d, a)) - 1
    uniroot(f, c(1e-9, 50))$root
  }, numeric(1))
  lines(rad * cos(th), rad * sin(th), col = cols[i], lwd = 2.6)
}
legend("topright", c("L1 (diamond)", expression(alpha == 0.25), expression(alpha == 0.5),
                     expression(alpha == 0.75), "L2 (disc)"),
       col = cols, lwd = 2.6, bty = "n", cex = 0.75)
mtext("corner on the axis = a residual pinned at exactly zero", side = 1,
      line = 3.4, cex = 0.62, col = "grey35")

# --- influence functions ----------------------------------------------------
u <- seq(-3, 3, length.out = 1600)
psi_q  <- sign(u)
psi_e  <- u
psi_el <- 2 * 0.5 * u + (1 - 0.5) * sign(u)          # elastile, alpha = 0.5, s = 1
cH <- 1
psi_H  <- pmax(pmin(u, cH), -cH)                     # Huber
psi_iv <- pmax(u, -cH)                               # one-sided inverted Huber
plot(NA, xlim = c(-3, 3), ylim = c(-3.2, 3.2),
     xlab = "residual u", ylab = expression(psi(u)),
     main = "Influence functions")
abline(h = 0, v = 0, col = "grey85")
# The inverted Huber coincides with Huber below the knot and with the expectile
# above it -- that IS its design, so draw the overlapping curves at different
# widths rather than pretending they are distinct.
lines(u, psi_e,  col = "#197a45", lwd = 4.5)
lines(u, psi_q,  col = "#1f5f8b", lwd = 4.5)
lines(u, psi_el, col = "#7d3c98", lwd = 2.6)
lines(u, psi_H,  col = "#a53a2b", lwd = 4.5, lty = 1)
lines(u, psi_iv, col = "#8a6d1f", lwd = 1.9, lty = 5)
legend("topleft", c("quantile (L1)", "expectile (L2)", expression(paste("elastile ", alpha == 0.5)),
                    "Huber", "inverted Huber (one-sided)"),
       col = c("#1f5f8b", "#197a45", "#7d3c98", "#a53a2b", "#8a6d1f"),
       lwd = c(4.5, 4.5, 2.6, 4.5, 1.9), lty = c(1, 1, 1, 1, 5), bty = "n", cex = 0.72)
text(1.9, -2.0, "inverted Huber = Huber below\nthe knot, expectile above",
     cex = 0.62, col = "#8a6d1f", adj = 0.5)
mtext("bounded psi = finite gross-error sensitivity; the mixture is unbounded for every alpha > 0",
      side = 1, line = 3.4, cex = 0.62, col = "grey35")
dev.off()
cat("wrote out/fig-loss-geometry.png\n")

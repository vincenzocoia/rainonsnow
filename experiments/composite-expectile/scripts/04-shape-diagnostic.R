# ---------------------------------------------------------------------------
# Why do the composite estimators return a negative shape at small n?
#
# In the n = 50 run every composite fit came back with xi < 0, though the truth
# is +0.2 and both estimators recover +0.2 from a huge sample. This traces the
# fitted shape as n grows, and measures how much of the weight sits at levels
# beyond the largest observation -- i.e. in pure extrapolation, where the
# empirical quantile/expectile function has saturated at the sample maximum and
# the loss can only pull the fitted tail down.
#
# Run with a wider shape bound than the main experiment so the bound itself
# cannot be the explanation.
#
# Output: out/shape-diagnostic.txt, out/fig-shape.png
# ---------------------------------------------------------------------------
source("R/setup.R")
source("R/config.R")
library(parallel)

XI_LO <<- -0.95        # widen: the main run bounds at -0.45
NS <- c(50, 100, 200, 500, 2000, 10000, 1e5)
NR <- 200

res <- list()
for (n in NS) {
  set.seed(777 + n)
  ds <- lapply(seq_len(NR), function(i) r_true(n))
  for (ty in c("quantile", "expectile")) {
    m <- do.call(rbind, mclapply(ds, function(y)
      fit_composite(y, GRID, W_FUN, ty, start = fit_lmom(y)), mc.cores = detectCores()))
    res[[sprintf("%s_%d", ty, n)]] <- m
  }
  cat("n =", n, "done\n")
}
saveRDS(res, "out/shape-diagnostic.rds")

# Share of the weight that lies above the largest observation, on average.
weight_beyond <- function(n) {
  tot <- sum(GRID$w_quad * W_FUN(GRID$p))
  set.seed(5); pmaxs <- replicate(400, p_true(max(r_true(n))))
  mean(sapply(pmaxs, function(pm) {
    k <- GRID$p > pm
    sum(GRID$w_quad[k] * W_FUN(GRID$p[k])) / tot
  }))
}

sink("out/shape-diagnostic.txt", split = TRUE)
cat("=== Fitted shape parameter vs sample size ===\n")
cat("Truth xi = 0.2. Shape bounded at", XI_LO, "(wider than the main run's -0.45).\n")
cat(NR, "replicates per cell.\n\n")
tab <- t(sapply(NS, function(n) {
  q <- sapply(c("quantile", "expectile"), function(ty) {
    x <- res[[sprintf("%s_%d", ty, n)]][, 3]
    c(median(x), mean(x <= XI_LO + 1e-6))
  })
  c(n = n, xi_L1 = q[1, 1], at_bound_L1 = q[2, 1],
    xi_L2 = q[1, 2], at_bound_L2 = q[2, 2],
    wt_beyond_data = weight_beyond(n))
}))
print(round(as.data.frame(tab), 4), row.names = FALSE)
cat("\nxi_*        : median fitted shape\n")
cat("at_bound_*  : fraction of fits sitting on the shape bound\n")
cat("wt_beyond_data: share of the weight function's mass at levels above the\n")
cat("                largest observation, where the empirical quantile and\n")
cat("                expectile functions have both saturated at the sample max\n")
cat("                and the loss can only pull the fitted tail downward.\n")
sink()

png("out/fig-shape.png", width = 1500, height = 620, res = 135)
par(mfrow = c(1, 2), mar = c(4.4, 4.6, 3.2, 1.2))
for (ty in c("quantile", "expectile")) {
  bx <- lapply(NS, function(n) res[[sprintf("%s_%d", ty, n)]][, 3])
  boxplot(bx, names = NS, outline = FALSE, log = "",
          xlab = "sample size n", ylab = expression(hat(xi)),
          main = sprintf("Composite %s (%s)", ty, if (ty == "quantile") "L1" else "L2"),
          col = adjustcolor(if (ty == "quantile") "#2980b9" else "#1e8449", 0.35))
  abline(h = 0.2, col = "black", lwd = 2, lty = 2)
  abline(h = XI_LO, col = "red", lty = 3)
  legend("bottomright", c("true xi = 0.2", "shape bound"), col = c("black", "red"),
         lty = c(2, 3), lwd = c(2, 1), bty = "n", cex = 0.8)
}
dev.off()
cat("\nwrote out/shape-diagnostic.txt and out/fig-shape.png\n")

# Copula transport, one conditional distribution per x.
#
# WHAT THIS SHOWS
#
# Each x gives its own estimate of the MARGINAL of Y, by transporting that x's
# conditional distribution back through the copula:
#
#     F_Y(y)  =  h^{-1}( F_{Y|X=x}(y) | u ),     u = F_X(x),
#
# with h(v|u) = dC(u,v)/du. Every x should give the same answer, so the spread
# across x is a direct, free diagnostic of whether the conditional model and the
# copula agree with each other.
#
# Three views, all as a function of x:
#   1. the whole transported marginal survival curve, one per x, against truth;
#   2. transported marginal quantiles at fixed exceedance probabilities;
#   3. the transported marginal EVI.
#
# For (3) the tail shape is fitted SEPARATELY at each x, so the flatness of the
# curve is a real diagnostic rather than an artefact. (Sharing one shape across x
# would make it flat by construction.) The transported marginal EVI is
# xi(x) / CEVI, and for the Gaussian copula CEVI = 1 - rho^2.
#
# DISTRIBUTIONAL LEARNING, NOT BINNING
#
# The conditional distributions come from a kernel-weighted conditional CDF: for
# a target x, every observation gets weight K((x_i - x)/h) and the conditional is
# the resulting weighted empirical distribution. That is the same kind of object
# a quantile regression forest returns -- a weighted empirical CDF over the
# training responses, forest weights instead of kernel weights -- so it feeds the
# same tail machinery the pipeline uses, and there is one distribution for every
# x rather than one per bin.
#
#   Rscript scripts/experiments/copula-transport-by-x.R
#
# Writes plots/copula-transport-by-x.png and prints a summary.

repo <- getwd()
while (!file.exists(file.path(repo, "DESCRIPTION")) && dirname(repo) != repo) {
  repo <- dirname(repo)
}
if (!file.exists(file.path(repo, "DESCRIPTION"))) {
  stop("Run this from inside the rainonsnow repository.", call. = FALSE)
}
source(file.path(repo, "R", "gpd_tail.R"))
source(file.path(repo, "R", "mixture_tail.R"))

set.seed(20240830)

N_OBS <- 2000L
RHO <- 0.7
XI_Y <- 0.25
BANDWIDTH <- 0.06 # kernel bandwidth on the u (rank) scale
TAIL_FRAC <- 0.35 # upper fraction of each conditional used for its GP tail
N_X <- 60L # x values at which to report
REPORT_P <- c(1e-2, 1e-3, 1e-4)

CEVI <- 1 - RHO^2
XI_COND <- CEVI * XI_Y

# ---------------------------------------------------------------------------
# Truth
# ---------------------------------------------------------------------------
SY <- function(y) (1 + XI_Y * y)^(-1 / XI_Y)
QY <- function(p) (p^(-XI_Y) - 1) / XI_Y

# Transport a conditional survival back to the marginal.
transport_survival <- function(s_cond, u, rho) {
  a <- stats::qnorm(s_cond, lower.tail = FALSE)
  stats::pnorm(a * sqrt(1 - rho^2) + rho * stats::qnorm(u), lower.tail = FALSE)
}

# ---------------------------------------------------------------------------
# Data
# ---------------------------------------------------------------------------
u <- stats::runif(N_OBS)
v <- stats::pnorm(RHO * stats::qnorm(u) + sqrt(1 - RHO^2) * stats::rnorm(N_OBS))
y <- QY(1 - v)
# x on an interpretable scale; the copula only sees ranks, so any monotone
# transform of u would do.
x <- stats::qnorm(u)

# ---------------------------------------------------------------------------
# Distributional learning: a weighted empirical conditional CDF at any x.
#
# Returns outcomes and weights -- exactly the shape of a quantregForest
# prediction -- so the tail extraction below is the same code the pipeline runs.
# ---------------------------------------------------------------------------
dl_conditional <- function(x0, x, y, h = BANDWIDTH) {
  # Work on the rank scale so the bandwidth means the same thing everywhere.
  ux <- stats::ecdf(x)(x)
  u0 <- stats::ecdf(x)(x0)
  w <- stats::dnorm((ux - u0) / h)
  keep <- w > 1e-8
  list(
    outcomes = y[keep],
    weights = w[keep] / sum(w[keep]),
    u = u0
  )
}

# Upper tail of a weighted empirical distribution, as in dl_tail_pieces().
weighted_tail <- function(cond, frac = TAIL_FRAC) {
  o <- order(cond$outcomes)
  ys <- cond$outcomes[o]
  ws <- cond$weights[o]
  cw <- cumsum(ws)
  thr <- ys[which(cw >= 1 - frac)[1]]
  above <- ys > thr
  if (sum(above) < 5L) {
    return(NULL)
  }
  wa <- ws[above]
  list(
    threshold = thr,
    excess = ys[above] - thr,
    weight = wa,
    tail_prob = sum(wa),
    n_eff = sum(wa)^2 / sum(wa^2)
  )
}

# ---------------------------------------------------------------------------
# One independent tail fit per x, then transport each to a marginal.
# ---------------------------------------------------------------------------
x_grid <- stats::quantile(x, seq(0.03, 0.97, length.out = N_X), names = FALSE)

fits <- lapply(x_grid, function(x0) {
  cond <- dl_conditional(x0, x, y)
  tl <- weighted_tail(cond)
  if (is.null(tl)) {
    return(NULL)
  }
  f <- fit_gpd_weighted(tl$excess, tl$weight)
  if (!is.finite(f$shape) || !is.finite(f$scale)) {
    return(NULL)
  }
  list(
    u = cond$u,
    threshold = tl$threshold,
    tail_prob = tl$tail_prob,
    scale = f$scale,
    shape = f$shape,
    n_eff = tl$n_eff
  )
})
ok <- !vapply(fits, is.null, logical(1L))
fits <- fits[ok]
x_ok <- x_grid[ok]

# Transported marginal survival implied by the conditional at one x.
marginal_survival_at_x <- function(fit, yy) {
  s_cond <- fit$tail_prob *
    gpd_survival(yy - fit$threshold, fit$scale, fit$shape)
  s_cond <- pmin(pmax(s_cond, 1e-300), 1 - 1e-12)
  transport_survival(s_cond, fit$u, RHO)
}

# Invert it for a quantile.
marginal_quantile_at_x <- function(fit, p) {
  vapply(
    p,
    function(pp) {
      lo <- fit$threshold
      hi <- lo + 1
      it <- 0L
      while (marginal_survival_at_x(fit, hi) > pp && it < 300L) {
        hi <- lo + (hi - lo) * 2
        it <- it + 1L
      }
      if (marginal_survival_at_x(fit, lo) < pp) {
        return(NA_real_)
      }
      for (i in seq_len(120L)) {
        mid <- (lo + hi) / 2
        if (marginal_survival_at_x(fit, mid) > pp) lo <- mid else hi <- mid
      }
      (lo + hi) / 2
    },
    numeric(1L)
  )
}

evi_by_x <- vapply(fits, function(f) f$shape, numeric(1L)) / CEVI
neff_by_x <- vapply(fits, function(f) f$n_eff, numeric(1L))
q_by_x <- vapply(
  fits,
  function(f) marginal_quantile_at_x(f, REPORT_P),
  numeric(length(REPORT_P))
)

# ---------------------------------------------------------------------------
# Combining across x.
#
# NOT by averaging the survival curves. The arithmetic mean of survivals is
# dominated by the heaviest x (the same failure mode as the unpooled mixture),
# and the geometric mean is dominated by the LIGHTEST, because log S runs to
# -infinity for a light tail. Both are badly biased; the geometric mean came out
# at 0.65x of the true 1e-4 return level here.
#
# Two defensible combinations instead:
#
#   median   -- the median across x of the per-x transported return level. Robust
#               to the wobble, makes no parametric assumption, but throws away
#               the fact that every x is estimating the same number.
#   shared   -- fit ONE conditional shape across all x (scales still free), then
#               transport. This is the estimator the pipeline should use: it pools
#               the tail information instead of averaging noisy separate answers.
# ---------------------------------------------------------------------------
conds <- lapply(x_grid, function(x0) {
  cond <- dl_conditional(x0, x, y)
  tl <- weighted_tail(cond)
  if (is.null(tl)) NULL else c(tl, list(u = cond$u))
})
conds <- conds[!vapply(conds, is.null, logical(1L))]

shared_fit <- fit_gpd_shared_shape(
  lapply(conds, `[[`, "excess"),
  lapply(conds, `[[`, "weight"),
  n_boot = 0L
)
shared_fits <- lapply(seq_along(conds), function(i) {
  list(
    u = conds[[i]]$u,
    threshold = conds[[i]]$threshold,
    tail_prob = conds[[i]]$tail_prob,
    scale = shared_fit$scale[i],
    shape = shared_fit$shape,
    n_eff = conds[[i]]$n_eff
  )
})
shared_fits <- shared_fits[vapply(shared_fits, function(f) is.finite(f$scale), logical(1L))]

q_shared <- vapply(
  shared_fits,
  function(f) marginal_quantile_at_x(f, REPORT_P),
  numeric(length(REPORT_P))
)

combined_median <- apply(q_by_x, 1, stats::median, na.rm = TRUE)
combined_shared <- apply(q_shared, 1, stats::median, na.rm = TRUE)

y_curve <- exp(seq(log(QY(0.3)), log(QY(1e-5)), length.out = 120L))

cat("\n=====================================================================\n")
cat("Copula transport with one conditional distribution per x\n")
cat("=====================================================================\n")
cat(sprintf(
  "n = %d, rho = %.2f, CEVI = %.2f\n",
  N_OBS,
  RHO,
  CEVI
))
cat(sprintf(
  "true marginal EVI %.3f, true conditional EVI %.3f\n",
  XI_Y,
  XI_COND
))
cat(sprintf("%d of %d x values gave a usable tail fit\n\n", length(fits), N_X))

cat(sprintf(
  "transported marginal EVI across x: median %.3f, IQR %.3f to %.3f\n",
  stats::median(evi_by_x),
  stats::quantile(evi_by_x, 0.25, names = FALSE),
  stats::quantile(evi_by_x, 0.75, names = FALSE)
))
inner <- x_ok > stats::quantile(x, 0.2) & x_ok < stats::quantile(x, 0.8)
cat(sprintf(
  "   middle 60%% of x: sd %.4f    outer 40%%: sd %.4f\n\n",
  stats::sd(evi_by_x[inner]),
  stats::sd(evi_by_x[!inner])
))

cat("marginal return level, truth vs combined estimate:\n")
for (j in seq_along(REPORT_P)) {
  truth_j <- QY(REPORT_P[j])
  spread <- stats::quantile(q_by_x[j, ], c(0.1, 0.9), na.rm = TRUE, names = FALSE)
  cat(sprintf(
    "   p = %-7.0e truth %8.2f   per-x 10-90%%: %7.2f to %7.2f\n",
    REPORT_P[j],
    truth_j,
    spread[1],
    spread[2]
  ))
}

cat(sprintf(
  "\nshared conditional shape across all x: %.4f (true %.4f)\n",
  shared_fit$shape,
  XI_COND
))
cat(sprintf(
  "  -> transported marginal EVI %.4f (true %.4f)\n\n",
  shared_fit$shape / CEVI,
  XI_Y
))
cat("combined marginal return levels:\n")
cat(sprintf(
  "%-10s %10s %10s %10s\n",
  "p",
  "truth",
  "median-x",
  "shared"
))
for (j in seq_along(REPORT_P)) {
  tj <- QY(REPORT_P[j])
  cat(sprintf(
    "%-10.0e %10.2f %8.2f%s %8.2f%s\n",
    REPORT_P[j],
    tj,
    combined_median[j],
    sprintf(" (%.2fx)", combined_median[j] / tj),
    combined_shared[j],
    sprintf(" (%.2fx)", combined_shared[j] / tj)
  ))
}

# ---------------------------------------------------------------------------
# Plots
# ---------------------------------------------------------------------------
plots_dir <- file.path(repo, "plots")
dir.create(plots_dir, showWarnings = FALSE, recursive = TRUE)
out_png <- file.path(plots_dir, "copula-transport-by-x.png")

grDevices::png(out_png, width = 1500, height = 1150, res = 140)
op <- graphics::par(
  mfrow = c(2, 2),
  mar = c(4.2, 4.4, 3.2, 1.2),
  cex.main = 1.0,
  col.axis = "grey30",
  col.lab = "grey25",
  fg = "grey45"
)

ramp <- grDevices::colorRampPalette(c("#2166ac", "#92c5de", "#f4a582", "#b2182b"))
cols <- ramp(length(fits))

# (1) one transported marginal survival curve per x
plot(
  NA,
  xlim = range(y_curve),
  ylim = c(1e-5, 0.5),
  log = "xy",
  xlab = "y",
  ylab = "P(Y > y)",
  main = "Transported marginal, one curve per x"
)
for (i in seq_along(fits)) {
  graphics::lines(
    y_curve,
    marginal_survival_at_x(fits[[i]], y_curve),
    col = grDevices::adjustcolor(cols[i], alpha.f = 0.55),
    lwd = 1
  )
}
graphics::lines(y_curve, SY(y_curve), col = "black", lwd = 3)
for (i in seq_along(shared_fits)) {
  graphics::lines(
    y_curve,
    marginal_survival_at_x(shared_fits[[i]], y_curve),
    col = grDevices::adjustcolor("#1a9850", alpha.f = 0.5),
    lwd = 1.2
  )
}
graphics::legend(
  "bottomleft",
  c("truth", "shared shape, per x", "own shape, per x"),
  col = c("black", "#1a9850", "#92c5de"),
  lwd = c(3, 1.2, 1),
  lty = c(1, 1, 1),
  bty = "n",
  cex = 0.85
)

# (2) transported marginal quantiles vs x
matplot_y <- t(q_by_x)
plot(
  NA,
  xlim = range(x_ok),
  ylim = range(matplot_y, na.rm = TRUE),
  log = "y",
  xlab = "x",
  ylab = "marginal return level",
  main = "Transported marginal quantiles vs x"
)
qcols <- c("#2166ac", "#c1663a", "#7b3294")
x_shared <- vapply(shared_fits, function(f) stats::qnorm(f$u), numeric(1L))
for (j in seq_along(REPORT_P)) {
  graphics::lines(
    x_ok,
    matplot_y[, j],
    col = grDevices::adjustcolor(qcols[j], alpha.f = 0.45),
    lwd = 1.6
  )
  graphics::lines(x_shared, q_shared[j, ], col = qcols[j], lwd = 2.6)
  graphics::abline(h = QY(REPORT_P[j]), col = "grey25", lty = 3, lwd = 1.2)
}
graphics::legend(
  "topleft",
  sprintf("p = %.0e", REPORT_P),
  col = qcols,
  lwd = 2.6,
  bty = "n",
  cex = 0.85
)
graphics::mtext(
  "bold = shared shape, faded = own shape per x, dotted = truth",
  side = 3,
  line = 0.2,
  cex = 0.66,
  col = "grey40"
)

# (3) transported marginal EVI vs x
plot(
  x_ok,
  evi_by_x,
  type = "l",
  lwd = 2,
  col = "#2166ac",
  xlab = "x",
  ylab = "transported marginal EVI",
  main = "Marginal EVI implied by each x",
  ylim = range(c(evi_by_x, XI_Y), na.rm = TRUE)
)
graphics::abline(h = XI_Y, col = "black", lwd = 2, lty = 2)
graphics::abline(h = shared_fit$shape / CEVI, col = "#1a9850", lwd = 2)
graphics::rug(
  x[x >= min(x_ok) & x <= max(x_ok)],
  col = grDevices::adjustcolor("grey40", alpha.f = 0.25)
)
graphics::legend(
  "topleft",
  c("own shape per x, / CEVI", "shared shape / CEVI", "true marginal EVI"),
  col = c("#2166ac", "#1a9850", "black"),
  lwd = 2,
  lty = c(1, 1, 2),
  bty = "n",
  cex = 0.8
)

# (4) how much information each x carries
plot(
  x_ok,
  neff_by_x,
  type = "h",
  col = "#8fc7d8",
  lwd = 3,
  xlab = "x",
  ylab = "effective sample size of the conditional tail",
  main = "Information behind each x"
)
graphics::lines(x_ok, neff_by_x, col = "#1f6f8b", lwd = 2)

graphics::par(op)
grDevices::dev.off()
cat(sprintf("\nWrote %s\n", out_png))
cat("=====================================================================\n")

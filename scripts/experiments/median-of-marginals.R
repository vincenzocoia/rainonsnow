# Combining per-x estimates of the marginal: why the pointwise MEDIAN.
#
# THE SETUP
#
# Each covariate value transports its own conditional tail back through the
# copula and so produces its own estimate of the WHOLE marginal of Y:
#
#     S_Y(y) = h^{-1}( S_{Y|X=x}(y) | u ),      u = F_X(x).
#
# These are repeated estimates of one quantity, not components of a mixture, so
# combining them is an averaging problem. That is what makes the median legal
# here and illegal for the operation it is often confused with: mixing the
# conditionals themselves is the law of total probability, and a median of
# survival curves is not an integral, so it has no business there. Once each
# curve is already an estimate of the same marginal, the median is exactly the
# right tool -- and a pointwise median of monotone survival curves is itself a
# monotone survival curve, so the answer is still a distribution.
#
# THE DATA-GENERATING PROCESS
#
# F_X and Y | X are chosen; the marginal of Y is derived from them (see
# R/transport_lab.R). Every conditional is exactly generalized Pareto. The
# marginal is not generalized Pareto at any level a sample reaches, so fitting
# one to Y directly is biased no matter how much data there is -- which is the
# comparison this experiment is against.
#
#   Rscript scripts/experiments/median-of-marginals.R
#
# Writes scripts/experiments/results/median-of-marginals-results.txt

repo <- getwd()
while (!file.exists(file.path(repo, "DESCRIPTION")) && dirname(repo) != repo) {
  repo <- dirname(repo)
}
if (!file.exists(file.path(repo, "DESCRIPTION"))) {
  stop("Run this from inside the rainonsnow repository.", call. = FALSE)
}
source(file.path(repo, "R", "gpd_tail.R"))
source(file.path(repo, "R", "copula_transport.R"))
source(file.path(repo, "R", "transport_lab.R"))

set.seed(20260829)

PS <- c(1e-3, 1e-4, 1e-5)
N_REP <- 60L
BW <- 0.03
N_ANCHOR <- 30L

dgp <- transport_dgp("pareto", alpha = 4, xi_cond = 0.1)
truth <- td_marginal_quantile(dgp, PS)
true_transport <- function(s, u) td_transport(dgp, s, u)

canned <- function(y) {
  g <- td_canned_gp(y)
  list(
    shape = g$shape,
    q = g$threshold + gpd_quantile_upper(PS / g$tail_prob, g$scale, g$shape)
  )
}

ratios <- function(m) paste(sprintf("%6.2f", m / truth), collapse = " ")

# ---------------------------------------------------------------------------
cat("\n=====================================================================\n")
cat("PART 1  The DGP, and why a canned fit cannot work on it\n")
cat("=====================================================================\n\n")
print(dgp)
cat("\nX ~ Pareto(4) is heavier than the conditional generalized Pareto, so\n")
cat("Breiman's lemma gives the marginal the covariate's index, 1/4, while every\n")
cat("conditional keeps 0.1. The ratio is the copula conditional EVI.\n\n")
qq <- td_marginal_quantile(dgp, 10^-(1:8))
cat(sprintf("%12s %12s %10s\n", "S_Y(y)", "y", "local EVI"))
for (i in seq_along(qq)) {
  cat(sprintf("%12.0e %12.2f %10.3f\n", 10^-i, qq[i], td_local_evi(dgp, qq[i])))
}
cat(sprintf(
  "\nA generalized Pareto fit lives in the top few percent of the sample, where\n"
))
cat(sprintf(
  "the local slope is %.3f to %.3f against an asymptote of %.3f. It reads that\n",
  td_local_evi(dgp, qq[2]),
  td_local_evi(dgp, qq[1]),
  td_asymptotic_evi(dgp)
))
cat("slope and then extrapolates it as if it had settled, so the error is in\n")
cat("the shape of the curve rather than in how precisely it was measured, and\n")
cat("more data does not remove it.\n")
utils::flush.console()

# ---------------------------------------------------------------------------
cat("\n\n=====================================================================\n")
cat("PART 2  One sample: the ensemble of per-x marginals\n")
cat("=====================================================================\n")

d <- td_simulate(dgp, 3000L)
cc <- td_conditionals(d$u, d$y, n_anchor = N_ANCHOR, bandwidth = BW)
me_shared <- transport_ensemble(cc$pieces, cc$u, family = true_transport)
me_free <- transport_ensemble(
  cc$pieces,
  cc$u,
  family = true_transport,
  shared_shape = FALSE
)

per_x <- function(me, p) {
  vapply(
    seq_along(me$u),
    function(j) {
      one <- me
      for (f in c("u", "threshold", "tail_prob", "scale", "shape", "weights")) {
        one[[f]] <- me[[f]][j]
      }
      marginal_ensemble_quantile(one, p)
    },
    numeric(1L)
  )
}

cat(sprintf("\nn = 3000, %d covariate values. True 1e-4 level: %.2f\n",
            length(me_shared$u), truth[2]))
cat("\nPer-x estimates of the 1e-4 level:\n")
cat(sprintf(
  "  one shape across x (%.3f):  min %.1f  median %.1f  max %.1f\n",
  me_shared$shape[1],
  min(per_x(me_shared, PS[2])),
  stats::median(per_x(me_shared, PS[2])),
  max(per_x(me_shared, PS[2]))
))
cat(sprintf(
  "  a shape per x (%.3f to %.3f): min %.1f  median %.1f  max %.1f\n",
  min(me_free$shape),
  max(me_free$shape),
  min(per_x(me_free, PS[2])),
  stats::median(per_x(me_free, PS[2])),
  max(per_x(me_free, PS[2]))
))
cat("\nThat spread is free diagnostic information: every x estimates the SAME\n")
cat("marginal, so disagreement between them measures how far the conditional\n")
cat("model and the copula are from being mutually consistent.\n")
utils::flush.console()

# ---------------------------------------------------------------------------
cat("\n\n=====================================================================\n")
cat("PART 3  The three combination rules, over repeated samples\n")
cat("        ratio of the estimate to the truth, at S = 1e-3, 1e-4, 1e-5\n")
cat("=====================================================================\n")

run <- function(n, shared) {
  out <- array(NA_real_, c(N_REP, 3L, 4L))
  for (r in seq_len(N_REP)) {
    s <- td_simulate(dgp, n)
    cc <- td_conditionals(s$u, s$y, n_anchor = N_ANCHOR, bandwidth = BW)
    me <- transport_ensemble(
      cc$pieces,
      cc$u,
      family = true_transport,
      shared_shape = shared
    )
    out[r, , 1] <- marginal_ensemble_quantile(me, PS, "median")
    out[r, , 2] <- marginal_ensemble_quantile(me, PS, "mean")
    out[r, , 3] <- marginal_ensemble_quantile(me, PS, "weighted_mean")
    out[r, , 4] <- canned(s$y)$q
  }
  out
}

for (n in c(1000L, 3000L, 10000L)) {
  for (shared in c(TRUE, FALSE)) {
    o <- run(n, shared)
    cat(sprintf(
      "\nn = %5d, %s\n",
      n,
      if (shared) "one shape across x" else "a shape per x"
    ))
    for (j in seq_len(4L)) {
      cat(sprintf(
        "   %-24s %s\n",
        c("median", "mean", "weighted mean", "canned GP on Y alone")[j],
        ratios(apply(o[, , j], 2L, stats::median, na.rm = TRUE))
      ))
    }
    utils::flush.console()
  }
}

cat("\nWith one shape across x the three rules are close, because the paths\n")
cat("differ only through their thresholds and scales. Let each x have its own\n")
cat("shape and they separate: the mean is dragged up by whichever x drew the\n")
cat("heaviest shape, which is the same domination that makes a mixture of\n")
cat("conditionals inherit max(xi_i), reappearing one level up. The median is\n")
cat("not, because a single heavy path moves it by one rank position.\n")
cat("\nThe canned fit degrades monotonically with how far it is extrapolated,\n")
cat("and more data does not rescue it -- it is fitting the wrong shape of curve,\n")
cat("not the right one imprecisely.\n")

# ---------------------------------------------------------------------------
cat("\n\n=====================================================================\n")
cat("PART 4  What decides whether the transport works\n")
cat("=====================================================================\n")

cat("\n4a. How hard the conditional shape error is amplified.\n")
cat("    Exact conditionals, true copula, shape shifted by delta. Nothing is\n")
cat("    estimated here, so what is left is the amplification alone.\n")

for (d in list(dgp, transport_dgp("lognormal", sdlog = 1, xi_cond = 0.2))) {
  tr <- td_marginal_quantile(d, PS)
  tp <- function(sc, u) td_transport(d, sc, u)
  u <- seq(0.03, 0.97, length.out = 30L)
  x <- td_x_quantile(d, u)
  span <- -log10(td_conditional_survival(d, tr[3], x))
  cat(sprintf(
    "\n    X ~ %s, CEVI = %.2f. For the 1e-5 event each conditional is pushed\n",
    if (d$x_dist == "pareto") sprintf("Pareto(%g)", d$alpha) else
      sprintf("LogNormal(0, %g)", d$sdlog),
    d$xi_cond / td_asymptotic_evi(d)
  ))
  cat(sprintf(
    "    %.1f to %.1f log units past its own body.\n",
    min(span),
    max(span)
  ))
  for (del in c(0, 0.005, 0.01, 0.02)) {
    me <- marginal_ensemble(
      u,
      rep(0, 30L),
      rep(1, 30L),
      x,
      d$xi_cond + del,
      family = tp
    )
    cat(sprintf(
      "      delta = %+.3f   %s\n",
      del,
      paste(sprintf("%6.2f", marginal_ensemble_quantile(me, PS) / tr),
            collapse = " ")
    ))
  }
}

cat("\n    delta = 0 returns exactly 1.00 everywhere, which is the identity the\n")
cat("    method rests on. Past that, the same shape error costs about half as\n")
cat("    much again where the span is wider. That span is what CEVI measures:\n")
cat("    at CEVI = 1 conditioning has not made the event any less extreme, so\n")
cat("    reaching a 1e-5 marginal event still needs a conditional event out at\n")
cat("    1e-13, and the estimate is y proportional to s^-xi the whole way.\n")

cat("\n\n4b. How well the conditional is learned, which matters more.\n")
cat("    The same ensemble, now with conditionals learned from the sample at\n")
cat("    four kernel bandwidths. n = 3000.\n")

for (d in list(dgp, transport_dgp("lognormal", sdlog = 1, xi_cond = 0.2))) {
  tr <- td_marginal_quantile(d, PS)
  tp <- function(sc, u) td_transport(d, sc, u)
  cat(sprintf(
    "\n    X ~ %s, CEVI = %.2f, true conditional shape %.3f\n",
    if (d$x_dist == "pareto") sprintf("Pareto(%g)", d$alpha) else
      sprintf("LogNormal(0, %g)", d$sdlog),
    d$xi_cond / td_asymptotic_evi(d),
    d$xi_cond
  ))
  for (bw in c(0.015, 0.03, 0.06, 0.12)) {
    set.seed(7)
    o <- t(replicate(N_REP, {
      s <- td_simulate(d, 3000L)
      cc <- td_conditionals(s$u, s$y, n_anchor = N_ANCHOR, bandwidth = bw)
      me <- transport_ensemble(cc$pieces, cc$u, family = tp)
      c(marginal_ensemble_quantile(me, PS, "median"), me$shape[1])
    }))
    cat(sprintf(
      "      bandwidth %.3f   shape %.3f   %s\n",
      bw,
      stats::median(o[, 4]),
      paste(sprintf("%6.2f",
                    apply(o[, 1:3], 2L, stats::median, na.rm = TRUE) / tr),
            collapse = " ")
    ))
    utils::flush.console()
  }
}

cat("\n    A kernel window on the rank scale spans a range of covariate values,\n")
cat("    so each learned conditional is itself a small scale mixture and comes\n")
cat("    out too heavy -- the same domination as everywhere else in this file,\n")
cat("    now inside a single conditional. How badly depends on how fast the\n")
cat("    conditional scale moves with the rank, which is much faster for the\n")
cat("    lognormal covariate: at bandwidth 0.12 its shape is 39% high against\n")
cat("    15% for the Pareto one, and the two effects compound.\n")
cat("\n    So the transport's accuracy is governed by the conditional model,\n")
cat("    and CEVI sets the price of getting it slightly wrong. Both regimes\n")
cat("    work at a bandwidth narrow enough to keep the conditional honest; at\n")
cat("    CEVI = 1 the margin for error is much thinner, and there is no tail\n")
cat("    heaviness for the decomposition to explain in the first place.\n")
cat("=====================================================================\n")

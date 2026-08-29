# Head to head: is decomposing across a predictor better than fitting the
# marginal tail directly?
#
# THE QUESTION THIS SETTLES
#
# That the transport beats MIXING the learned conditionals is not in doubt --
# the mixture carries the conditional tail index and is structurally too light.
# The open question is whether it beats the obvious alternative: ignore the
# covariate and fit one generalized Pareto to Y. This is the comparison the
# whole decomposition argument stands on, so the competitor is given every
# advantage that can be given without making it clairvoyant, and one that does:
# an oracle that picks the threshold in hindsight, per replicate and per target.
#
# THE SWEEP
#
# The covariate's informativeness is varied while the marginal's asymptotic tail
# index is held fixed at 1/4. With X ~ Pareto(4) and Y | X = x generalized
# Pareto with scale x and shape xi, Breiman's lemma gives the marginal the
# covariate's index whenever xi < 1/4, so
#
#     marginal EVI = 0.25 throughout,        CEVI = xi / 0.25.
#
# Sweeping xi therefore sweeps CEVI from 0.1 (the covariate explains nearly all
# the tail heaviness) to 0.8 (it explains little), against a fixed target. The
# derived marginal is recomputed by quadrature at each setting; nothing about it
# is assumed. xi = 0.25 is left out: there E[W^4] diverges, Breiman no longer
# applies, and the index picks up a logarithmic factor.
#
# WHAT IS ESTIMATED, AND WHAT IS NOT
#
# The kernel bandwidth follows one fixed rule, 0.1487 * n^(-1/5), calibrated
# once to give 0.03 at n = 3000 and never tuned per setting -- otherwise the
# transport would be scored on a choice the competitor is not allowed to make.
# The copula is fitted from the sample's Kendall tau in two of the three
# transport rows; the third uses the process's own copula and is an upper bound,
# not a method, isolating how much of the error is conditional estimation rather
# than copula misspecification.
#
#   Rscript scripts/experiments/transport-vs-marginal-fit.R
#
# Writes scripts/experiments/results/transport-vs-marginal-fit-results.txt

repo <- getwd()
while (!file.exists(file.path(repo, "DESCRIPTION")) && dirname(repo) != repo) {
  repo <- dirname(repo)
}
if (!file.exists(file.path(repo, "DESCRIPTION"))) {
  stop("Run this from inside the rainonsnow repository.", call. = FALSE)
}
source(file.path(repo, "R", "gpd_tail.R"))
source(file.path(repo, "R", "mixture_tail.R"))
source(file.path(repo, "R", "copula_transport.R"))
source(file.path(repo, "R", "transport_lab.R"))

ALPHA <- 4
XI_GRID <- c(0.025, 0.05, 0.10, 0.15, 0.20) # CEVI 0.1 to 0.8
N_GRID <- c(1000L, 3000L, 10000L)
N_REP <- 60L
N_ANCHOR <- 30L
PS <- c(1e-4, 1e-5)
GP_THRESHOLDS <- c(0.80, 0.90, 0.95, 0.975)

bandwidth_for <- function(n) 0.1487 * n^(-1 / 5)

METHODS <- c(
  "transport, true copula (upper bound)",
  "transport, fitted Gaussian copula",
  "transport, fitted survival Clayton",
  "mixture over covariate values",
  "GP on Y alone, 90% threshold",
  "GP on Y alone, oracle threshold"
)

# One replicate: every estimator on the same sample, so the comparison is paired.
one_rep <- function(dgp, n, truth) {
  s <- td_simulate(dgp, n)
  cc <- td_conditionals(
    s$u,
    s$y,
    n_anchor = N_ANCHOR,
    bandwidth = bandwidth_for(n)
  )
  tau <- stats::cor(s$u, s$y, method = "kendall")
  out <- matrix(NA_real_, length(METHODS), length(PS))

  me <- try(
    transport_ensemble(
      cc$pieces,
      cc$u,
      family = function(a, b) td_transport(dgp, a, b)
    ),
    silent = TRUE
  )
  if (!inherits(me, "try-error")) {
    out[1, ] <- marginal_ensemble_quantile(me, PS, "median")

    g <- me
    g$family <- "gaussian"
    g$par <- sin(pi / 2 * tau)
    out[2, ] <- marginal_ensemble_quantile(g, PS, "median")

    sc <- me
    sc$family <- "survival_clayton"
    sc$par <- if (tau > 0) 2 * tau / (1 - tau) else NA_real_
    if (is.finite(sc$par)) {
      out[3, ] <- marginal_ensemble_quantile(sc, PS, "median")
    }

    # The law of total probability over the covariate values, which is what the
    # pipeline currently does. Equal weights are right because the anchors are
    # evenly spaced in u.
    mt <- try(
      mixture_tail(me$threshold, me$tail_prob, me$scale, me$shape),
      silent = TRUE
    )
    if (!inherits(mt, "try-error")) {
      out[4, ] <- mixture_tail_quantile(mt, PS)
    }
  }

  # The competitor, at the conventional threshold and at the best one available
  # in hindsight -- chosen per replicate and per target probability, which no
  # real analysis could do.
  cand <- vapply(
    GP_THRESHOLDS,
    function(p0) {
      g <- td_canned_gp(s$y, p = p0)
      g$threshold + gpd_quantile_upper(PS / g$tail_prob, g$scale, g$shape)
    },
    numeric(length(PS))
  )
  out[5, ] <- cand[, GP_THRESHOLDS == 0.90]
  out[6, ] <- vapply(
    seq_along(PS),
    function(i) {
      r <- abs(log(cand[i, ] / truth[i]))
      if (all(!is.finite(r))) NA_real_ else cand[i, which.min(r)]
    },
    numeric(1L)
  )
  out
}

summarise <- function(est, truth) {
  # Median ratio for bias, root mean square log ratio for total error, and the
  # share of replicates that produced nothing at all.
  r <- est / truth
  ok <- is.finite(r) & r > 0
  list(
    med = if (any(ok)) stats::median(r[ok]) else NA_real_,
    rmse = if (any(ok)) sqrt(mean(log(r[ok])^2)) else NA_real_,
    na = mean(!ok)
  )
}

results <- array(
  NA_real_,
  c(length(METHODS), length(PS), length(XI_GRID), length(N_GRID), 3L)
)

cat("\n=====================================================================\n")
cat("Copula transport against a generalized Pareto fitted to Y alone\n")
cat("=====================================================================\n")
cat(sprintf(
  "\nX ~ Pareto(%g); Y | X = x ~ GP(scale = x, shape = xi). The marginal is\n",
  ALPHA
))
cat("derived by quadrature and has asymptotic tail index 0.25 throughout, so\n")
cat("sweeping xi sweeps CEVI against a fixed target. Kernel bandwidth follows\n")
cat("one fixed rule and is never tuned. Paired samples, ")
cat(sprintf("%d replicates each.\n", N_REP))

for (ni in seq_along(N_GRID)) {
  n <- N_GRID[ni]
  for (xi in seq_along(XI_GRID)) {
    dgp <- transport_dgp("pareto", alpha = ALPHA, xi_cond = XI_GRID[xi])
    truth <- td_marginal_quantile(dgp, PS)
    set.seed(1000L * ni + xi)
    est <- array(NA_real_, c(N_REP, length(METHODS), length(PS)))
    for (r in seq_len(N_REP)) {
      est[r, , ] <- one_rep(dgp, n, truth)
    }
    for (m in seq_along(METHODS)) {
      for (p in seq_along(PS)) {
        sm <- summarise(est[, m, p], truth[p])
        results[m, p, xi, ni, ] <- c(sm$med, sm$rmse, sm$na)
      }
    }
    cat(sprintf(
      "  ... n = %5d, CEVI = %.1f done\n",
      n,
      XI_GRID[xi] / 0.25
    ))
    utils::flush.console()
  }
}

# ---------------------------------------------------------------------------
block <- function(ni, p, stat, label) {
  cat(sprintf("\n%s, n = %d, target P(Y > y) = %s\n", label, N_GRID[ni],
              format(PS[p], scientific = TRUE)))
  cat(sprintf("  %-38s", "CEVI"))
  cat(paste(sprintf("%8.1f", XI_GRID / 0.25), collapse = ""), "\n")
  for (m in seq_along(METHODS)) {
    cat(sprintf("  %-38s", METHODS[m]))
    cat(paste(sprintf("%8.2f", results[m, p, , ni, stat]), collapse = ""), "\n")
  }
}

cat("\n\n=====================================================================\n")
cat("MEDIAN RATIO TO THE TRUE RETURN LEVEL   (1.00 is unbiased)\n")
cat("=====================================================================\n")
for (ni in seq_along(N_GRID)) {
  for (p in seq_along(PS)) {
    block(ni, p, 1L, "median ratio")
  }
}

cat("\n\n=====================================================================\n")
cat("ROOT MEAN SQUARE LOG RATIO   (total error; lower is better)\n")
cat("=====================================================================\n")
for (ni in seq_along(N_GRID)) {
  for (p in seq_along(PS)) {
    block(ni, p, 2L, "rms log ratio")
  }
}

cat("\n\n=====================================================================\n")
cat("SHARE OF REPLICATES RETURNING NOTHING\n")
cat("=====================================================================\n")
cat("\nA transport row fails when the fitted conditional shape is negative\n")
cat("enough that its transported survival never reaches the target: the\n")
cat("estimator declines to extrapolate rather than guessing. That is a\n")
cat("property worth seeing, not a defect to hide, but it also means the rows\n")
cat("above are summarised over the replicates that did return something.\n")
for (ni in seq_along(N_GRID)) {
  block(ni, length(PS), 3L, "failure rate")
}

saveRDS(
  list(
    results = results,
    methods = METHODS,
    cevi = XI_GRID / 0.25,
    n = N_GRID,
    probs = PS
  ),
  file.path(repo, "scripts", "experiments", "results",
            "transport-vs-marginal-fit.rds")
)
cat("\n=====================================================================\n")

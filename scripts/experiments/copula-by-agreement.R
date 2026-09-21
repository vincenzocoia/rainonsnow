# Choosing the copula by making the covariate values agree.
#
# WHY
#
# The head-to-head (transport-vs-marginal-fit.R) found that the transport is
# exact given the true copula and the worst estimator in the table given one
# fitted from the sample -- worse than the mixture it was meant to replace, and
# worse than ignoring the covariate. The whole gap is copula model risk, it does
# not shrink with n, and it gets worse as CEVI falls, because the more of the
# tail the covariate could explain the more of the answer comes from the corner
# of the copula. Kendall's tau is a statistic of the BODY. The corner is what
# matters and the body barely constrains it.
#
# But the transport supplies a criterion the marginal fit does not have. Every
# covariate value transports to an estimate of the SAME marginal, so under the
# right copula they must agree. Disagreement is an observable, computable
# without the truth, and it is a function of the copula. That suggests choosing
# the copula to minimise it -- an estimator that targets the corner directly
# rather than inheriting it from a body statistic.
#
# This script asks two questions.
#
#   A. Does the disagreement DETECT the misspecification? If it does, the method
#      self-diagnoses, which no marginal fit can do, and that is worth something
#      even if the point estimate is not.
#
#   B. Does minimising it REPAIR the estimate?
#
# THE DISAGREEMENT STATISTIC
#
# The standard deviation of log S across covariate values, summed over three
# reference levels taken from the sample itself (the 99th and 99.9th percentiles
# and the maximum). No truth enters it, so it is available in practice.
#
#   Rscript scripts/experiments/copula-by-agreement.R
#
# Writes scripts/experiments/results/copula-by-agreement-results.txt

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

ALPHA <- 4
XI_GRID <- c(0.025, 0.05, 0.10, 0.15, 0.20) # CEVI 0.1 to 0.8
N_GRID <- c(3000L, 10000L)
N_REP <- 60L
N_ANCHOR <- 30L
PS <- c(1e-4, 1e-5)

bandwidth_for <- function(n) 0.1487 * n^(-1 / 5)

# The candidate copulas the criterion chooses among.
CANDIDATES <- rbind(
  data.frame(family = "gaussian", par = seq(0.05, 0.95, by = 0.05)),
  data.frame(
    family = "clayton",
    par = exp(seq(log(0.05), log(30), length.out = 19L))
  ),
  data.frame(
    family = "survival_clayton",
    par = exp(seq(log(0.05), log(30), length.out = 19L))
  )
)

swap <- function(me, family, par) {
  me$family <- family
  me$par <- par
  me
}

# Between-covariate disagreement, computed without the truth.
disagreement <- function(me, y_ref) {
  paths <- marginal_ensemble_paths(me, y_ref)
  v <- apply(paths, 1L, function(s) {
    s <- s[is.finite(s) & s > 0]
    if (length(s) < 3L) NA_real_ else stats::sd(log(s))
  })
  if (all(is.na(v))) NA_real_ else sum(v, na.rm = TRUE)
}

one_rep <- function(dgp, n, truth) {
  s <- td_simulate(dgp, n)
  cc <- td_conditionals(
    s$u,
    s$y,
    n_anchor = N_ANCHOR,
    bandwidth = bandwidth_for(n)
  )
  tp <- function(a, b) td_transport(dgp, a, b)
  me <- try(transport_ensemble(cc$pieces, cc$u, family = tp), silent = TRUE)
  if (inherits(me, "try-error")) {
    return(NULL)
  }
  y_ref <- c(stats::quantile(s$y, c(0.99, 0.999), names = FALSE), max(s$y))
  y_ref <- pmax(y_ref, marginal_ensemble_lower(me) * 1.001)

  tau <- stats::cor(s$u, s$y, method = "kendall")
  me_tau <- swap(me, "gaussian", sin(pi / 2 * tau))

  # The criterion, over every candidate.
  d <- vapply(
    seq_len(nrow(CANDIDATES)),
    function(i) {
      disagreement(swap(me, CANDIDATES$family[i], CANDIDATES$par[i]), y_ref)
    },
    numeric(1L)
  )
  best <- which.min(d)
  me_best <- swap(me, CANDIDATES$family[best], CANDIDATES$par[best])

  g <- td_canned_gp(s$y)
  list(
    est = rbind(
      `transport, true copula (upper bound)` =
        marginal_ensemble_quantile(me, PS, "median"),
      `transport, copula from Kendall tau` =
        marginal_ensemble_quantile(me_tau, PS, "median"),
      `transport, copula by minimum disagreement` =
        marginal_ensemble_quantile(me_best, PS, "median"),
      `GP on Y alone, 90% threshold` =
        g$threshold + gpd_quantile_upper(PS / g$tail_prob, g$scale, g$shape)
    ),
    disagree = c(
      true = disagreement(me, y_ref),
      tau = disagreement(me_tau, y_ref),
      best = d[best]
    ),
    chosen = as.character(CANDIDATES$family[best]),
    chosen_par = CANDIDATES$par[best]
  )
}

cat("\n=====================================================================\n")
cat("Choosing the copula by making the covariate values agree\n")
cat("=====================================================================\n")
cat("\nSame process as transport-vs-marginal-fit.R: X ~ Pareto(4),\n")
cat("Y | X = x generalized Pareto with scale x, marginal index 0.25 throughout,\n")
cat("CEVI swept by the conditional shape. The disagreement statistic uses only\n")
cat(sprintf("the sample. %d replicates per cell.\n", N_REP))

for (n in N_GRID) {
  bias <- rmse <- matrix(NA_real_, 4L, length(XI_GRID))
  dis <- matrix(NA_real_, 3L, length(XI_GRID))
  corr <- numeric(length(XI_GRID))
  picked <- character(length(XI_GRID))

  for (xi in seq_along(XI_GRID)) {
    dgp <- transport_dgp("pareto", alpha = ALPHA, xi_cond = XI_GRID[xi])
    truth <- td_marginal_quantile(dgp, PS)
    set.seed(4000L + 100L * which(N_GRID == n) + xi)
    est <- array(NA_real_, c(N_REP, 4L, length(PS)))
    dd <- matrix(NA_real_, N_REP, 3L)
    ch <- character(N_REP)
    for (r in seq_len(N_REP)) {
      o <- one_rep(dgp, n, truth)
      if (is.null(o)) {
        next
      }
      est[r, , ] <- o$est
      dd[r, ] <- o$disagree
      ch[r] <- o$chosen
    }
    # Report at the deeper target; the shallower one is in the saved object.
    p <- length(PS)
    for (m in 1:4) {
      r <- est[, m, p] / truth[p]
      ok <- is.finite(r) & r > 0
      bias[m, xi] <- stats::median(r[ok])
      rmse[m, xi] <- sqrt(mean(log(r[ok])^2))
    }
    dis[, xi] <- apply(dd, 2L, stats::median, na.rm = TRUE)
    # Does the disagreement predict the error, replicate by replicate, for the
    # copula an analyst would actually have fitted?
    e <- abs(log(est[, 2L, p] / truth[p]))
    ok <- is.finite(e) & is.finite(dd[, 2L]) & dd[, 2L] > 0
    corr[xi] <- stats::cor(log(dd[ok, 2L]), e[ok], method = "spearman")
    picked[xi] <- names(sort(table(ch[nzchar(ch)]), decreasing = TRUE))[1]
    cat(sprintf("  ... n = %5d, CEVI = %.1f done\n", n, XI_GRID[xi] / 0.25))
    utils::flush.console()
  }

  hdr <- function(lab) {
    cat(sprintf("\n%s, n = %d, target P(Y > y) = %s\n", lab, n,
                format(PS[length(PS)], scientific = TRUE)))
    cat(sprintf("  %-44s", "CEVI"))
    cat(paste(sprintf("%8.1f", XI_GRID / 0.25), collapse = ""), "\n")
  }
  nm <- c(
    "transport, true copula (upper bound)",
    "transport, copula from Kendall tau",
    "transport, copula by minimum disagreement",
    "GP on Y alone, 90% threshold"
  )

  hdr("MEDIAN RATIO TO TRUTH")
  for (m in 1:4) {
    cat(sprintf("  %-44s", nm[m]))
    cat(paste(sprintf("%8.2f", bias[m, ]), collapse = ""), "\n")
  }

  hdr("ROOT MEAN SQUARE LOG RATIO")
  for (m in 1:4) {
    cat(sprintf("  %-44s", nm[m]))
    cat(paste(sprintf("%8.2f", rmse[m, ]), collapse = ""), "\n")
  }

  hdr("DISAGREEMENT BETWEEN COVARIATE VALUES (median; 0 is perfect)")
  for (m in 1:3) {
    cat(sprintf("  %-44s",
                c("under the true copula", "under the Kendall tau copula",
                  "under the chosen copula")[m]))
    cat(paste(sprintf("%8.2f", dis[m, ]), collapse = ""), "\n")
  }

  cat(sprintf("\n  %-44s", "Spearman corr(log disagreement, |log error|)"))
  cat(paste(sprintf("%8.2f", corr), collapse = ""), "\n")
  cat(sprintf("  %-44s", "family most often chosen"))
  cat(paste(sprintf("%8s", substr(picked, 1, 7)), collapse = ""), "\n")
  utils::flush.console()
}
cat("\n=====================================================================\n")

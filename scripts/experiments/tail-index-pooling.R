# Experiment: how should the tail index of a cell's runoff marginal be estimated?
#
# THE PROBLEM
#
# The marginal runoff distribution for a cell is the equal-weight mixture of its
# peak-hour predictive distributions, each grafted to its own GP tail. If
# component i has shape xi_i then its survival is regularly varying with index
# 1/xi_i, so the mixture's index is min_i(1/xi_i): the cell's tail is set by
# max_i(xi_i) alone. Each xi_i is fitted to the upper half of one predictive
# distribution, whose effective sample size is small, so max_i(xi_i) is the
# maximum of a few hundred noisy estimates and is inflated far above any of the
# underlying truth.
#
# WHAT THIS SCRIPT DOES
#
# Part 1 measures the damage and compares five estimators within a single cell.
# Part 2 asks whether shapes should also be smoothed across neighbouring cells.
#
# Runs on base R alone -- no probaverse, no tidyverse -- so it can be checked
# quickly and in isolation from the pipeline. A few minutes at the defaults.
#
#   Rscript scripts/experiments/tail-index-pooling.R

# Walk up from the working directory to the package root, so this runs from
# anywhere in the repo (the rest of the pipeline expects the root).
repo <- getwd()
while (!file.exists(file.path(repo, "DESCRIPTION")) && dirname(repo) != repo) {
  repo <- dirname(repo)
}
if (!file.exists(file.path(repo, "DESCRIPTION"))) {
  stop("Run this from inside the rainonsnow repository.", call. = FALSE)
}
source(file.path(repo, "R", "gpd_tail.R"))
source(file.path(repo, "R", "mixture_tail.R"))
source(file.path(repo, "R", "tail_shape_field.R"))

set.seed(20240828)

N_REP <- 40L # replicate cells per configuration
K_HOURS <- 150L # peak hours per cell
M_TAIL <- 15L # effective tail sample per predictive distribution
XI_TRUE <- 0.15 # true shape, shared within a cell
EVENTS_PER_YEAR <- 3
# Spanning three decades of return period, because the way the mixture inherits
# max_i(xi_i) is an asymptotic statement: it barely shows at T = 10 and
# dominates by T = 1000.
REPORT_T <- c(10, 100, 1000)

# ---------------------------------------------------------------------------
# Simulate one cell.
#
# Every peak hour has its own predictive distribution: same tail shape, a scale
# that depends on that hour's rain/snowmelt. That is exactly the assumption a
# shared shape encodes, and it is the standard non-stationary POT model --
# covariates move the scale, not the index.
# ---------------------------------------------------------------------------
simulate_cell <- function(k = K_HOURS, m = M_TAIL, xi = XI_TRUE) {
  scale <- stats::runif(k, 1, 6)
  threshold <- stats::runif(k, 8, 12)
  excess <- lapply(scale, function(s) {
    gpd_quantile_upper(stats::runif(m), s, xi)
  })
  list(
    scale = scale,
    threshold = threshold,
    tail_prob = rep(0.5, k),
    excess = excess
  )
}

# Truth for the cell: the mixture built from the true parameters.
true_mixture <- function(cell, xi = XI_TRUE) {
  mixture_tail(cell$threshold, cell$tail_prob, cell$scale, xi)
}

# ---------------------------------------------------------------------------
# The estimators.
#
# Each returns a per-component (scale, shape) pair, from which the mixture and
# its return levels follow. The differences are entirely in how the shape is
# obtained.
# ---------------------------------------------------------------------------
estimators <- list(
  # Status quo: fit each predictive distribution's tail on its own.
  independent = function(cell) {
    fits <- lapply(cell$excess, fit_gpd_weighted)
    list(
      scale = vapply(fits, `[[`, numeric(1L), "scale"),
      shape = vapply(fits, `[[`, numeric(1L), "shape")
    )
  },

  # Fit independently, then force the median shape on everyone and refit the
  # scales. Cheap, robust, no joint optimisation -- a reasonable fallback.
  median_shape = function(cell) {
    fits <- lapply(cell$excess, fit_gpd_weighted)
    xi <- stats::median(vapply(fits, `[[`, numeric(1L), "shape"), na.rm = TRUE)
    scale <- vapply(
      cell$excess,
      function(y) {
        gpd_profile_scale(y, rep(1 / length(y), length(y)), xi)$scale
      },
      numeric(1L)
    )
    list(scale = scale, shape = rep(xi, length(scale)))
  },

  # Pooled maximum likelihood: one shape, per-component scales.
  pooled = function(cell) {
    fit <- fit_gpd_shared_shape(cell$excess, n_boot = 0L)
    list(scale = fit$scale, shape = rep(fit$shape, length(fit$scale)))
  },

  # Pooled, with the incidental-parameter bias removed by bootstrap.
  pooled_bc = function(cell) {
    fit <- fit_gpd_shared_shape(cell$excess, n_boot = 25L)
    list(scale = fit$scale, shape = rep(fit$shape, length(fit$scale)))
  },

  # Reference point: the shape handed to us exactly right, scales still fitted.
  # Nothing can beat this, so it bounds how much the shape matters at all.
  oracle_shape = function(cell) {
    scale <- vapply(
      cell$excess,
      function(y) {
        gpd_profile_scale(y, rep(1 / length(y), length(y)), XI_TRUE)$scale
      },
      numeric(1L)
    )
    list(scale = scale, shape = rep(XI_TRUE, length(scale)))
  }
)

# ---------------------------------------------------------------------------
# Part 1: within one cell.
# ---------------------------------------------------------------------------
flush_cat <- function(...) { cat(...); utils::flush.console() }

cat("\n=====================================================================\n")
cat("PART 1  Estimating the shape within a cell\n")
cat(sprintf(
  "        %d peak hours, %d effective tail points each, true xi = %.2f\n",
  K_HOURS,
  M_TAIL,
  XI_TRUE
))
cat(sprintf("        %d replicate cells\n", N_REP))
cat("=====================================================================\n\n")

target_p <- 1 / (REPORT_T * EVENTS_PER_YEAR)

results <- list()
for (nm in names(estimators)) {
  shape_hat <- numeric(N_REP)
  rl <- matrix(NA_real_, N_REP, length(REPORT_T))
  rl_true <- matrix(NA_real_, N_REP, length(REPORT_T))

  for (r in seq_len(N_REP)) {
    cell <- simulate_cell()
    est <- estimators[[nm]](cell)

    # The mixture's effective index is set by the heaviest component.
    shape_hat[r] <- max(est$shape, na.rm = TRUE)

    mt <- mixture_tail(cell$threshold, cell$tail_prob, est$scale, est$shape)
    rl[r, ] <- mixture_tail_quantile(mt, target_p)
    rl_true[r, ] <- mixture_tail_quantile(true_mixture(cell), target_p)
  }

  flush_cat(sprintf("   %-14s done\n", nm))
  results[[nm]] <- list(
    shape = shape_hat,
    ratio = rl / rl_true
  )
}

cat(sprintf(
  "%-14s %17s   %s\n",
  "estimator",
  "effective xi",
  "return level / truth (median, sd)"
))
cat(sprintf(
  "%-14s %8s %8s   %s\n",
  "",
  "mean",
  "bias",
  paste(sprintf("%17s", paste0("T=", REPORT_T, "y")), collapse = "")
))
cat(strrep("-", 80), "\n")
for (nm in names(results)) {
  rr <- results[[nm]]
  cat(sprintf(
    "%-14s %8.3f %+8.3f   %s\n",
    nm,
    mean(rr$shape, na.rm = TRUE),
    mean(rr$shape, na.rm = TRUE) - XI_TRUE,
    paste(
      sprintf(
        "%9.2fx %5.2f",
        apply(rr$ratio, 2, stats::median, na.rm = TRUE),
        apply(rr$ratio, 2, stats::sd, na.rm = TRUE)
      ),
      collapse = ""
    )
  ))
}

cat("\nReading this table:\n")
cat("  'effective xi' is max_i(xi_i), the index the mixture actually inherits.\n")
cat("  'return level / truth' compares against the mixture built from the TRUE\n")
cat("  per-hour parameters, so 1.00x is exactly right.\n")

# The return-level column mixes two independent error sources, and reporting
# only the total hides which is which:
#
#   * shape error  -- the mixture inherits max_i(xi_i).
#   * scale noise  -- each hour's scale comes from a handful of points, and the
#                     return level is convex in the scales, so noise alone
#                     pushes the mixture up (Jensen). `oracle_shape` is given
#                     the true shape and still misses, which measures this term
#                     on its own.
oracle_infl <- apply(results$oracle_shape$ratio, 2, stats::median, na.rm = TRUE)
cat("\nDecomposition:\n")
cat(sprintf(
  "  Scale noise alone (oracle_shape, correct shape given): %s\n",
  paste(sprintf("T=%dy %.2fx", REPORT_T, oracle_infl), collapse = "   ")
))
cat("  Any estimator's remaining gap from those numbers is its shape error.\n")
cat("\n  Note how `independent` looks respectable at T=10 and falls apart as T\n")
cat("  grows. That is not accuracy: its typical shape is biased low, which\n")
cat("  offsets the scale inflation, while max_i(xi_i) takes over further out.\n")
cat("  Two large errors cancelling at one return period is not a method.\n")

# ---------------------------------------------------------------------------
# Part 2: should shapes be smoothed across neighbouring cells?
#
# The shape is now one number per cell, but that number is still noisy, and the
# true shape varies smoothly over terrain. This asks whether borrowing from
# neighbours beats treating each cell on its own -- and whether it is safe when
# the field varies more sharply than assumed.
# ---------------------------------------------------------------------------
cat("\n\n=====================================================================\n")
cat("PART 2  Smoothing the per-cell shape across neighbours\n")
cat("=====================================================================\n\n")

GRID <- 5L # GRID x GRID cells
N_FIELD <- 6L # replicate fields

field_experiment <- function(gradient, label) {
  g <- expand.grid(x = seq_len(GRID), y = seq_len(GRID))
  # A smooth underlying shape field: a linear trend across the grid.
  truth <- XI_TRUE + gradient * (g$x - mean(g$x))

  err_raw <- err_smooth <- err_global <- numeric(N_FIELD)
  for (r in seq_len(N_FIELD)) {
    est <- se <- numeric(nrow(g))
    for (i in seq_len(nrow(g))) {
      cell <- simulate_cell(k = 40L, xi = truth[i])
      fit <- fit_gpd_shared_shape(cell$excess, n_boot = 10L)
      est[i] <- fit$shape
      se[i] <- fit$shape_se
    }
    sm <- smooth_tail_shape(est, se, g$x, g$y, radius = 1)
    err_raw[r] <- mean((est - truth)^2)
    err_smooth[r] <- mean((sm$shape - truth)^2)
    err_global[r] <- mean((mean(est) - truth)^2)
  }

  flush_cat(sprintf("%s (shape ranges %.2f to %.2f across the grid)\n",
    label, min(truth), max(truth)))
  cat(sprintf(
    "   RMSE  per-cell %.4f   neighbour-smoothed %.4f   single global %.4f\n",
    sqrt(mean(err_raw)),
    sqrt(mean(err_smooth)),
    sqrt(mean(err_global))
  ))
  cat(sprintf(
    "   smoothing changes RMSE by %+.1f%%, pooling everything by %+.1f%%\n\n",
    100 * (sqrt(mean(err_smooth)) / sqrt(mean(err_raw)) - 1),
    100 * (sqrt(mean(err_global)) / sqrt(mean(err_raw)) - 1)
  ))
}

field_experiment(0.00, "Flat field   -- shape genuinely constant")
field_experiment(0.02, "Gentle trend -- shape varies slowly")
field_experiment(0.06, "Steep trend  -- shape varies sharply")

cat("=====================================================================\n")
cat("Part 2 checks the risk as well as the benefit: if the true field is\n")
cat("steep, over-smoothing would show up as smoothing losing to per-cell\n")
cat("estimates, and pooling to a single global shape losing badly.\n")
cat("=====================================================================\n")

# Does decomposing across a driver reduce bias at HYDROLOGICAL sample sizes?
#
# WHAT CHANGED FROM THE EARLIER EXPERIMENTS
#
# They ran at n = 1000 to 10000 and targeted exceedance probabilities of 1e-4 to
# 1e-5. Hydrology has 50 to 100 years of record and cares about the 20- to
# 200-year event. Both of those choices were wrong for the question, and the
# second one especially: 1e-5 is a 100,000-year event, which is where
# extrapolation error compounds and where no method deserves to be believed.
# Here the record IS the sample size, the peaks are annual maxima, and the
# T-year level is the 1/T upper quantile.
#
# THE PROCESS
#
# See R/two_process_dgp.R. An annual peak is the larger of a snowmelt-driven
# component, bounded above because a snowpack is finite, and a rain-driven one
# whose scale is set by the observed rainfall and which is heavy-tailed. The
# marginal has a kink at the bounded component's endpoint and belongs to no
# standard family; the conditionals do not, and the covariate selects which
# regime the peak is in. That is the decomposition argument at its strongest,
# and it is why the conditional SHAPE has to move with the covariate here, not
# just the scale.
#
# THE ESTIMATORS
#
# Standard practice is a distribution fitted to the peak itself. Against it:
# the drivers fitted separately and recombined (which needs them observed); the
# two-driver structure fitted from (rainfall, peak) pairs alone and then either
# integrated over a fitted rainfall distribution or averaged over the observed
# rainfall values; and the same fitted conditionals transported through a
# copula, with the copula estimated from Kendall's tau and from the upper part
# of the sample only.
#
# The tau-versus-upper contrast is the point of Coia, Joe and Nolde (2024): the
# conditional tail index is governed by the copula near the top, and the
# restriction here is on v alone across all of u, not to the upper-right corner,
# because for families with a non-constant conditional extreme value index the
# relevant behaviour is spread along that edge.
#
#   Rscript scripts/experiments/hydrological-decomposition.R
#
# Writes scripts/experiments/results/hydrological-decomposition-results.txt

repo <- getwd()
while (!file.exists(file.path(repo, "DESCRIPTION")) && dirname(repo) != repo) {
  repo <- dirname(repo)
}
if (!file.exists(file.path(repo, "DESCRIPTION"))) {
  stop("Run this from inside the rainonsnow repository.", call. = FALSE)
}
for (f in c("gpd_tail", "copula_transport", "two_process_dgp",
            "hydro_estimators")) {
  source(file.path(repo, "R", paste0(f, ".R")))
}

TT <- c(20, 50, 100, 200)
PS <- 1 / TT
N_GRID <- c(50L, 100L, 200L, 500L, 2000L)
N_REP <- 200L
POT_FRAC <- 0.5 # top half, which is all n = 50 can spare

dgp <- two_process_dgp()
TRUTH <- tp_marginal_quantile(dgp, PS)

METHODS <- c(
  "GEV on the peak (standard practice)",
  "GP-POT on the peak, top 50%",
  "POT regression on rainfall, shape varies",
  "POT regression on rainfall, one shared shape",
  "structure fitted, integrated over fitted F_X",
  "structure fitted, averaged over observed x",
  "drivers observed, each fitted marginally",
  "drivers observed, rain driver fitted on rainfall",
  "transport, true copula",
  "transport, Gaussian from Kendall tau",
  "transport, Gaussian from the upper half",
  "transport, survival Clayton, upper half"
)
# Assign by name so adding a row cannot silently shift another one.
ROW <- stats::setNames(seq_along(METHODS), METHODS)

# A monotone interpolator for the TRUE conditional at one x, so the true-copula
# transport does not need a bisection per evaluation.
true_cond_inverse <- function(x) {
  ly <- seq(base::log(1e-3), base::log(1e7), length.out = 600L)
  ls <- base::log(tp_conditional_survival(dgp, exp(ly), x))
  ok <- is.finite(ls) & c(TRUE, diff(ls) < 0)
  function(s) {
    exp(stats::approx(rev(ls[ok]), rev(ly[ok]), xout = base::log(s),
                      rule = 2)$y)
  }
}

# Median across observations of each one's transported marginal survival.
transport_curve <- function(s_cond_fun, u, family, par, y) {
  vapply(
    y,
    function(yy) {
      s <- s_cond_fun(yy)
      stats::median(
        cop_h_surv_inverse(pmin(pmax(s, 1e-300), 1 - 1e-12), u, family, par),
        na.rm = TRUE
      )
    },
    numeric(1L)
  )
}

one_rep <- function(n) {
  d <- tp_simulate(dgp, n)
  out <- matrix(NA_real_, length(METHODS), length(PS))

  # 1. GEV on the peak.
  out[ROW[["GEV on the peak (standard practice)"]], ] <-
    gev_return_level(fit_gev(d$y), PS)

  # 2. Generalized Pareto above the top half of the peaks.
  thr <- stats::quantile(d$y, 1 - POT_FRAC, names = FALSE)
  e <- d$y[d$y > thr] - thr
  if (length(e) >= 5L) {
    g <- fit_gpd_weighted(e, rep(1, length(e)))
    if (is.finite(g$shape)) {
      out[ROW[["GP-POT on the peak, top 50%"]], ] <-
        thr + gpd_quantile_upper(PS / POT_FRAC, g$scale, g$shape)
    }
  }

  # 3. Both drivers observed, each fitted on its own, recombined by
  #    F_Y = F_A F_B. The decomposition in its original form.
  fa <- fit_gp_ml(d$a)
  fb <- fit_gp_ml(d$b)
  if (is.finite(fa$shape) && is.finite(fb$shape)) {
    out[ROW[["drivers observed, each fitted marginally"]], ] <- invert_survival(
      function(y) {
        sa <- gpd_survival(y, fa$scale, fa$shape)
        sb <- gpd_survival(y, fb$scale, fb$shape)
        sa + sb - sa * sb
      },
      PS
    )
  }

  # A generic covariate model: no knowledge of how the peak is assembled, just a
  # non-stationary peaks-over-threshold regression on rainfall, removed over the
  # fitted rainfall distribution. The two variants are the shared-shape question
  # this repository currently answers with gp_tail.shape_pooling.
  fxr <- fit_gp_ml(d$x)
  for (vs in c(TRUE, FALSE)) {
    pr <- fit_pot_regression(d$x, d$y, fxr, varying_shape = vs)
    if (isTRUE(pr$converged)) {
      out[ROW[[if (vs) {
        "POT regression on rainfall, shape varies"
      } else {
        "POT regression on rainfall, one shared shape"
      }]], ] <- invert_survival(function(y) potr_marginal_survival(pr, y), PS)
    }
  }

  # Both drivers observed AND correctly specified: A on its own, the rain
  #    driver on the rainfall that scaled it. The upper bound for a
  #    decomposition, and the thing row 3 falls short of.
  fbx <- fit_scaled_gp(d$x, d$b)
  fx0 <- fit_gp_ml(d$x)
  if (is.finite(fa$shape) && is.finite(fbx$shape) && is.finite(fx0$shape)) {
    ftrue_form <- list(
      scale_a = fa$scale, shape_a = fa$shape,
      b_coef = fbx$b_coef, shape_b = fbx$shape
    )
    out[ROW[["drivers observed, rain driver fitted on rainfall"]], ] <-
      invert_survival(function(y) tpf_marginal_survival(ftrue_form, fx0, y), PS)
  }

  # The structure fitted from (x, y) alone. Everything below shares it.
  f <- fit_two_process(d$x, d$y)
  if (!is.finite(f$shape_b)) {
    return(out)
  }
  fx <- fit_gp_ml(d$x)

  # 4/5. Integrate over the fitted rainfall distribution, or average over the
  #      observed rainfall values.
  out[ROW[["structure fitted, integrated over fitted F_X"]], ] <-
    invert_survival(function(y) tpf_marginal_survival(f, fx, y), PS)
  out[ROW[["structure fitted, averaged over observed x"]], ] <-
    invert_survival(function(y) tpf_mixture_survival(f, d$x, y), PS)

  # 6-9. The same conditionals, transported instead of integrated.
  u <- rank(d$x) / (n + 1)
  v <- rank(d$y) / (n + 1)
  s_cond <- function(yy) tpf_conditional_survival(f, yy, d$x)
  ygrid <- exp(seq(base::log(TRUTH[1] * 0.2), base::log(TRUTH[4] * 20),
                   length.out = 90L))

  # This row builds one interpolator of the true conditional per observation, so
  # it is capped at 200 of them; a median over 200 draws is plenty.
  idx <- if (n > 200L) sample.int(n, 200L) else seq_len(n)
  inv <- lapply(d$x[idx], true_cond_inverse)
  true_curve <- vapply(
    ygrid,
    function(yy) {
      sc <- pmin(pmax(tpf_conditional_survival(f, yy, d$x[idx]), 1e-300),
                 1 - 1e-12)
      yy_true <- vapply(seq_along(idx), function(i) inv[[i]](sc[i]), numeric(1L))
      stats::median(tp_marginal_survival(dgp, yy_true), na.rm = TRUE)
    },
    numeric(1L)
  )
  out[ROW[["transport, true copula"]], ] <- interp_return(ygrid, true_curve, PS)

  tau <- stats::cor(u, v, method = "kendall")
  rho_tau <- sin(pi / 2 * tau)
  out[ROW[["transport, Gaussian from Kendall tau"]], ] <- interp_return(
    ygrid, transport_curve(s_cond, u, "gaussian", rho_tau, ygrid), PS
  )
  rho_up <- fit_copula_upper(u, v, "gaussian", 0.5)
  if (is.finite(rho_up)) {
    out[ROW[["transport, Gaussian from the upper half"]], ] <- interp_return(
      ygrid, transport_curve(s_cond, u, "gaussian", rho_up, ygrid), PS
    )
  }
  th_up <- fit_copula_upper(u, v, "survival_clayton", 0.5)
  if (is.finite(th_up)) {
    out[ROW[["transport, survival Clayton, upper half"]], ] <- interp_return(
      ygrid, transport_curve(s_cond, u, "survival_clayton", th_up, ygrid), PS
    )
  }
  out
}

# Read return levels off a survival curve evaluated on a grid.
interp_return <- function(y, s, p) {
  ok <- is.finite(s) & s > 0 & c(TRUE, diff(s) < 0)
  if (sum(ok) < 5L) {
    return(rep(NA_real_, length(p)))
  }
  out <- exp(stats::approx(base::log(s[ok]), base::log(y[ok]),
                           xout = base::log(p), rule = 1)$y)
  out
}

cat("\n=====================================================================\n")
cat("Decomposition at hydrological sample sizes and return periods\n")
cat("=====================================================================\n\n")
print(dgp)
cat(sprintf(
  "\nTrue return levels, T = 20 / 50 / 100 / 200 years: %s\n",
  paste(sprintf("%.2f", TRUTH), collapse = "  ")
))
cat(sprintf(
  "Local tail index there: %s (asymptotically %.2f)\n",
  paste(sprintf("%.3f", tp_local_evi(dgp, TRUTH)), collapse = " "),
  max(dgp$shape_x, dgp$shape_b)
))
cat(sprintf("\n%d replicates per sample size.\n", N_REP))

store <- array(NA_real_, c(length(METHODS), length(PS), length(N_GRID), 3L))
for (ni in seq_along(N_GRID)) {
  n <- N_GRID[ni]
  est <- array(NA_real_, c(N_REP, length(METHODS), length(PS)))
  set.seed(7000L + ni)
  for (r in seq_len(N_REP)) {
    est[r, , ] <- one_rep(n)
  }
  for (m in seq_along(METHODS)) {
    for (p in seq_along(PS)) {
      ratio <- est[, m, p] / TRUTH[p]
      ok <- is.finite(ratio) & ratio > 0
      store[m, p, ni, 1] <- if (any(ok)) stats::median(ratio[ok]) else NA_real_
      store[m, p, ni, 2] <- if (any(ok)) {
        sqrt(mean(base::log(ratio[ok])^2))
      } else {
        NA_real_
      }
      store[m, p, ni, 3] <- mean(!ok)
    }
  }
  cat(sprintf("  ... n = %5d done\n", n))
  utils::flush.console()
}

show <- function(stat, label, fmt) {
  cat("\n\n=====================================================================\n")
  cat(label, "\n")
  cat("=====================================================================\n")
  for (ni in seq_along(N_GRID)) {
    cat(sprintf("\nn = %d years\n", N_GRID[ni]))
    cat(sprintf("  %-48s", "T (years)"))
    cat(paste(sprintf("%8d", TT), collapse = ""), "\n")
    for (m in seq_along(METHODS)) {
      cat(sprintf("  %-48s", METHODS[m]))
      cat(paste(sprintf(fmt, store[m, , ni, stat]), collapse = ""), "\n")
    }
  }
}

show(1L, "MEDIAN RATIO TO THE TRUE RETURN LEVEL   (1.00 is unbiased)", "%8.2f")
show(2L, "ROOT MEAN SQUARE LOG RATIO   (total error; lower is better)", "%8.2f")
show(3L, "SHARE OF REPLICATES RETURNING NOTHING", "%8.2f")

saveRDS(
  list(store = store, methods = METHODS, T = TT, n = N_GRID, truth = TRUTH),
  file.path(repo, "scripts", "experiments", "results",
            "hydrological-decomposition.rds")
)
cat("\n=====================================================================\n")

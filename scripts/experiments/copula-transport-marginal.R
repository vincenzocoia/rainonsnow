# Experiment: does going through the copula give a better marginal tail for Y
# than estimating that tail directly?
#
# BACKGROUND
#
# Write U = F_X(X), V = F_Y(Y), and let W = 1/(1-V), which is standard Pareto
# (EVI 1) marginally. The *copula conditional EVI*, CEVI(u), is the extreme value
# index of W given U = u. It is a property of the copula alone, and for heavy or
# light tailed Y it satisfies
#
#     EVI(Y | X = x)  =  CEVI(u) * EVI(Y),        u = F_X(x),
#
# with CEVI in [0, 1]. So every conditional distribution has a LIGHTER tail than
# the marginal, and the copula says by exactly how much. For the Gaussian copula
# CEVI = 1 - rho^2, constant in u (verified in part 1 below): the stronger the
# dependence, the more of the marginal's tail heaviness is explained by the
# predictor rather than by the conditional noise.
#
# TWO CONSEQUENCES, WHICH THIS SCRIPT TESTS
#
# 1. The pipeline's marginal is the equal-weight mixture of the fitted
#    conditionals over the OBSERVED covariates. Every component shares the
#    conditional EVI, so the finite mixture is regularly varying with that index
#    -- which is CEVI * EVI(Y), strictly lighter than the truth whenever
#    CEVI < 1. The true marginal only gets its full heaviness from covariate
#    values beyond any finite sample. So the mixture marginal is asymptotically
#    too light by construction, no matter how well each conditional is fitted.
#
# 2. Transporting a conditional back through the copula,
#       F_Y(y) = h^{-1}( F_{Y|X=x}(y) | u ),      h(v|u) = dC(u,v)/du,
#    recovers the marginal exactly, and splits the estimation into two pieces
#    that are each better determined than EVI(Y) is directly: a conditional tail
#    index (from all n points, at a shallower extrapolation) and a copula
#    parameter (from all n ranks, not just the tail).
#
# THE COMPARISON
#
# The mixture and the transport are given IDENTICAL conditional fits, so the only
# thing that differs is how those fits are aggregated into a marginal. Direct POT
# on the y values alone is the external baseline.
#
# Base R only. Runs in a few minutes.
#
#   Rscript scripts/experiments/copula-transport-marginal.R

repo <- getwd()
while (!file.exists(file.path(repo, "DESCRIPTION")) && dirname(repo) != repo) {
  repo <- dirname(repo)
}
if (!file.exists(file.path(repo, "DESCRIPTION"))) {
  stop("Run this from inside the rainonsnow repository.", call. = FALSE)
}
source(file.path(repo, "R", "gpd_tail.R"))
source(file.path(repo, "R", "mixture_tail.R"))

set.seed(20240829)

# Sourcing this file with SWEEP_MODE set defines the estimators without running
# the report, so scripts/experiments/copula-transport-sweep.R can reuse them.
SWEEP_MODE <- exists("SWEEP_MODE") && isTRUE(SWEEP_MODE)

N_REP <- 300L # replicate datasets
N_OBS <- 300L # observations per dataset (POT peaks in one cell)
RHO <- 0.7 # Gaussian copula parameter
XI_Y <- 0.25 # true marginal EVI
N_BINS <- 8L # covariate bins for the conditional fits
COND_TAIL_FRAC <- 0.4 # fraction of each bin used for its GP tail
POT_FRAC <- 0.2 # fraction of the sample used by direct POT
REPORT_P <- c(1e-2, 1e-3, 1e-4) # marginal exceedance probabilities

flush_cat <- function(...) {
  cat(...)
  utils::flush.console()
}

# ---------------------------------------------------------------------------
# Copula pieces (Gaussian), all in log space where the tail needs it.
# ---------------------------------------------------------------------------
h_fun <- function(v, u, rho) {
  stats::pnorm((stats::qnorm(v) - rho * stats::qnorm(u)) / sqrt(1 - rho^2))
}

# Survival version of the transport: given S_cond = P(Y > y | U = u), return
# P(Y > y). Working from survival probabilities directly keeps the far tail
# accurate where 1 - F would round to 1.
transport_survival <- function(s_cond, u, rho) {
  a <- stats::qnorm(s_cond, lower.tail = FALSE) # = qnorm(1 - s_cond)
  stats::pnorm(
    a * sqrt(1 - rho^2) + rho * stats::qnorm(u),
    lower.tail = FALSE
  )
}

# True marginal and conditional for the copula DGP; used by part 2 and part 3.
SY <- function(y) (1 + XI_Y * y)^(-1 / XI_Y)
QY <- function(p) (p^(-XI_Y) - 1) / XI_Y
S_cond_true <- function(y, u) 1 - h_fun(1 - SY(y), u, RHO)

# ---------------------------------------------------------------------------
# Part 1: CEVI of the Gaussian copula is 1 - rho^2.
# ---------------------------------------------------------------------------
if (!SWEEP_MODE) {
cat("\n=====================================================================\n")
cat("PART 1  The copula conditional EVI (CEVI) of the Gaussian copula\n")
cat("=====================================================================\n\n")

log_surv_W <- function(log_w, u, rho) {
  z <- stats::qnorm(-log_w, lower.tail = FALSE, log.p = TRUE)
  stats::pnorm(
    (z - rho * stats::qnorm(u)) / sqrt(1 - rho^2),
    lower.tail = FALSE,
    log.p = TRUE
  )
}
cat(sprintf("%8s %8s %12s %12s\n", "rho", "u", "CEVI", "1 - rho^2"))
for (rho in c(0.3, 0.5, 0.7, 0.9)) {
  for (u in c(0.2, 0.5, 0.9)) {
    lw <- log(10) * c(60, 80)
    slope <- diff(log_surv_W(lw, u, rho)) / diff(lw)
    cat(sprintf("%8.1f %8.1f %12.5f %12.5f\n", rho, u, -1 / slope, 1 - rho^2))
  }
}
cat("\nCEVI is constant in u for this copula, and every conditional tail is\n")
cat("lighter than the marginal by exactly this factor.\n")

# ---------------------------------------------------------------------------
# Part 2: what tail does the finite mixture over observed covariates have?
# ---------------------------------------------------------------------------
cat("\n\n=====================================================================\n")
cat("PART 2  The tail of the mixture over OBSERVED covariates\n")
cat("=====================================================================\n\n")

p_probe <- 10^-seq(1, 8, length.out = 8)
y_probe <- QY(p_probe)
local_evi <- function(s, y) -1 / (diff(log(s)) / diff(log(y)))

cat(sprintf(
  "true marginal EVI %.3f, true conditional EVI %.3f (= CEVI x marginal)\n\n",
  XI_Y,
  (1 - RHO^2) * XI_Y
))
cat(sprintf("%-28s %s\n", "local EVI at p =", paste(
  sprintf("%7.0e", p_probe[-1]),
  collapse = " "
)))
cat(sprintf("%-28s %s\n", "true marginal", paste(
  sprintf("%7.3f", local_evi(SY(y_probe), y_probe)),
  collapse = " "
)))
for (n in c(200L, 2000L, 20000L)) {
  u <- (seq_len(n) - 0.5) / n
  s_mix <- vapply(y_probe, function(yy) mean(S_cond_true(yy, u)), numeric(1L))
  cat(sprintf("%-28s %s\n", sprintf("mixture over n = %d", n), paste(
    sprintf("%7.3f", local_evi(s_mix, y_probe)),
    collapse = " "
  )))
}
cat("\nThe mixture tracks the marginal for a while and then peels away toward\n")
cat("the conditional EVI. Where it peels away depends on n, and at the sample\n")
cat("sizes this project has, that happens inside the range of interest.\n")
}

# ---------------------------------------------------------------------------
# Part 3: estimator comparison on simulated data.
# ---------------------------------------------------------------------------
if (!SWEEP_MODE) {
  cat("\n\n=====================================================================\n")
  cat("PART 3  Estimating the marginal return level from n observations\n")
  cat("=====================================================================\n\n")
}

# --- data generating processes ---------------------------------------------
# A: Gaussian copula + GPD marginal. CEVI = 1 - rho^2 < 1, so the conditionals
#    are genuinely lighter than the marginal. This is the regime the transport
#    is designed for.
dgp_copula <- function(n) {
  u <- stats::runif(n)
  v <- stats::pnorm(
    RHO * stats::qnorm(u) + sqrt(1 - RHO^2) * stats::rnorm(n)
  )
  list(u = u, y = QY(1 - v), true_S = SY)
}

# B: conditional GPD with a bounded covariate effect on the scale, common shape.
#    Here the marginal EVI EQUALS the conditional EVI (CEVI = 1), so transport
#    has nothing to gain and this checks that it does not do harm.
dgp_scale <- function(n) {
  u <- stats::runif(n)
  sigma <- exp(1.2 * (u - 0.5)) # bounded in u, so no extra tail heaviness
  y <- gpd_quantile_upper(stats::runif(n), sigma, XI_Y)
  ug <- (seq_len(2001L) - 0.5) / 2001L
  sg <- exp(1.2 * (ug - 0.5))
  true_S <- function(yy) {
    vapply(yy, function(t) mean(gpd_survival(t, sg, XI_Y)), numeric(1L))
  }
  list(u = u, y = y, true_S = true_S)
}

# C: Gaussian copula again, but the marginal of Y is a mixture of two Pareto
#    tails with different indices. It is still regularly varying with EVI 0.25,
#    but the lighter component keeps contributing well into the tail, so the
#    marginal is NOT close to a GP at any level a sample of this size reaches.
#
#    This is the case the whole argument is about. DGP A hands direct POT an
#    exactly-GP marginal, which is the best case it could possibly be given and
#    is not what a mixture over covariates produces. Here the marginal has the
#    "nuance" that a single GP fit to the top few observations has to bulldoze
#    through, while each conditional remains well behaved.
XI_HEAVY <- 0.25
XI_LIGHT <- 0.05
S_nuanced <- function(y) {
  0.5 * (1 + y)^(-1 / XI_HEAVY) + 0.5 * (1 + y)^(-1 / XI_LIGHT)
}
# Quantile by interpolation on a dense log grid: needed thousands of times per
# replicate, so root-finding each call would dominate the runtime.
.nuance_grid_y <- c(0, exp(seq(log(1e-6), log(1e12), length.out = 6000L)))
.nuance_grid_s <- S_nuanced(.nuance_grid_y)
Q_nuanced <- function(p) {
  stats::approx(
    x = log(.nuance_grid_s),
    y = log1p(.nuance_grid_y),
    xout = log(p),
    rule = 2
  )$y |>
    expm1()
}

dgp_nuanced <- function(n) {
  u <- stats::runif(n)
  v <- stats::pnorm(
    RHO * stats::qnorm(u) + sqrt(1 - RHO^2) * stats::rnorm(n)
  )
  list(u = u, y = Q_nuanced(1 - v), true_S = S_nuanced)
}

# --- shared conditional fitting step ----------------------------------------
# Bin the covariate, then fit each bin's upper tail with ONE shape shared across
# bins. Both the mixture and the transport start from exactly this.
# Two ways to parameterise the per-bin scales:
#
#   "free"      -- one scale per bin. This is what the pipeline does per peak
#                  hour, and it is a Neyman-Scott setup: the number of nuisance
#                  scales grows with the number of bins while each is fitted from
#                  a handful of points, so the shared shape comes out biased low.
#   "loglinear" -- log sigma is linear in u, so the whole conditional family has
#                  three parameters instead of B + 1. This is the same fix as
#                  making sigma a smooth function of the covariates in the real
#                  pipeline, and it matters much more here because the transport
#                  divides the shape by CEVI and so amplifies any error in it.
fit_conditionals <- function(
  u,
  y,
  n_bins = N_BINS,
  frac = COND_TAIL_FRAC,
  scale_model = c("free", "loglinear"),
  bias_correct = FALSE
) {
  scale_model <- match.arg(scale_model)
  br <- stats::quantile(u, seq(0, 1, length.out = n_bins + 1L), names = FALSE)
  br[1] <- -Inf
  br[length(br)] <- Inf
  bin <- cut(u, br, labels = FALSE)

  thr <- num <- ubar <- numeric(n_bins)
  excess <- weight <- vector("list", n_bins)
  for (b in seq_len(n_bins)) {
    yb <- y[bin == b]
    ubar[b] <- mean(u[bin == b])
    num[b] <- length(yb)
    if (length(yb) < 5L) {
      excess[[b]] <- numeric(0)
      weight[[b]] <- numeric(0)
      next
    }
    thr[b] <- stats::quantile(yb, 1 - frac, names = FALSE)
    e <- yb[yb > thr[b]] - thr[b]
    excess[[b]] <- e
    weight[[b]] <- rep(1, length(e))
  }

  if (scale_model == "free") {
    fit <- fit_gpd_shared_shape(
      excess,
      weight,
      n_boot = if (bias_correct) 20L else 0L,
      bias_correct = bias_correct
    )
    shape <- fit$shape
    scale <- fit$scale
  } else {
    fit <- fit_gpd_loglinear_scale(excess, ubar)
    shape <- fit$shape
    scale <- fit$scale
  }

  list(
    threshold = thr,
    scale = scale,
    shape = shape,
    tail_prob = vapply(excess, length, integer(1L)) / pmax(num, 1L),
    u = ubar,
    n = num
  )
}

# Shared shape with log sigma_b = a + b * u_b, profiled over the shape.
fit_gpd_loglinear_scale <- function(excess, u, shape_bounds = c(-0.45, 1)) {
  keep <- vapply(excess, length, integer(1L)) >= 2L
  ex <- excess[keep]
  uu <- u[keep]
  out <- list(shape = NA_real_, scale = rep(NA_real_, length(excess)))
  if (length(ex) < 2L) {
    return(out)
  }

  nll <- function(par, xi) {
    sg <- exp(par[1] + par[2] * uu)
    total <- 0
    for (i in seq_along(ex)) {
      v <- gpd_nllh(ex[[i]], rep(1, length(ex[[i]])), sg[i], xi)
      if (!is.finite(v)) {
        return(Inf)
      }
      total <- total + v
    }
    total
  }

  # The start must satisfy the support constraint sigma > -xi * max(excess) when
  # xi < 0, or optim() is handed an Inf at the initial point and cannot move.
  mx <- max(unlist(ex))
  mn <- mean(unlist(ex))
  start_for <- function(xi) {
    c(log(max(mn, if (xi < 0) -xi * mx * 1.5 else 0)), 0)
  }
  profile <- function(xi) {
    stats::optim(start_for(xi), nll, xi = xi, control = list(reltol = 1e-8))$value
  }
  opt <- stats::optimize(profile, shape_bounds, tol = 1e-4)
  fin <- stats::optim(
    start_for(opt$minimum),
    nll,
    xi = opt$minimum,
    control = list(reltol = 1e-8)
  )
  out$shape <- opt$minimum
  out$scale[keep] <- exp(fin$par[1] + fin$par[2] * uu)
  out
}

# A local (k-nearest-neighbour) conditional model, as a stand-in for the
# quantile regression forest the real pipeline uses.
#
# Binning is a poor proxy for a fitted conditional model, and the bin width turns
# out to drive the answer: wide bins make each "conditional" a mixture over the
# bin's covariate range, which is heavier than the conditional at its centre and
# transports through as an inflated marginal; narrow bins starve each fit. A kNN
# neighbourhood is centred on its anchor and the neighbourhoods overlap, so the
# same data supports many more conditional fits. That is much closer to what a
# forest does, and it is the version whose result should be believed.
fit_conditionals_knn <- function(
  u,
  y,
  n_anchor = 20L,
  k = max(60L, floor(length(y) / 8)),
  frac = COND_TAIL_FRAC
) {
  anchors <- seq(0.5 / n_anchor, 1 - 0.5 / n_anchor, length.out = n_anchor)
  thr <- num <- numeric(n_anchor)
  excess <- weight <- vector("list", n_anchor)
  k <- min(k, length(y))

  for (a in seq_len(n_anchor)) {
    nb <- order(abs(u - anchors[a]))[seq_len(k)]
    yb <- y[nb]
    num[a] <- k
    thr[a] <- stats::quantile(yb, 1 - frac, names = FALSE)
    e <- yb[yb > thr[a]] - thr[a]
    excess[[a]] <- e
    # Neighbourhoods overlap, so each fit reuses data. Down-weighting by the
    # overlap factor keeps the pooled likelihood from counting each observation
    # n_anchor * k / n times over.
    weight[[a]] <- rep(length(y) / (n_anchor * k), length(e))
  }

  fit <- fit_gpd_shared_shape(excess, weight, n_boot = 0L)
  list(
    threshold = thr,
    scale = fit$scale,
    shape = fit$shape,
    tail_prob = vapply(excess, length, integer(1L)) / pmax(num, 1L),
    u = anchors,
    n = num
  )
}

# --- the three estimators ---------------------------------------------------

# Direct POT on y alone: the standard approach, ignoring the covariate.
est_direct_pot <- function(u, y, p) {
  k <- max(10L, floor(POT_FRAC * length(y)))
  thr <- stats::quantile(y, 1 - k / length(y), names = FALSE)
  e <- y[y > thr] - thr
  f <- fit_gpd_weighted(e, rep(1, length(e)))
  if (!is.finite(f$shape)) {
    return(rep(NA_real_, length(p)))
  }
  zeta <- length(e) / length(y)
  thr + gpd_quantile_upper(p / zeta, f$scale, f$shape)
}

# The pipeline: equal-weight mixture of the fitted conditionals.
est_dl_mixture <- function(cf, p) {
  ok <- is.finite(cf$threshold) &
    is.finite(cf$scale) &
    cf$tail_prob > 0 &
    cf$n > 0
  if (!any(ok) || !is.finite(cf$shape)) {
    return(rep(NA_real_, length(p)))
  }
  mt <- mixture_tail(
    threshold = cf$threshold[ok],
    tail_prob = cf$tail_prob[ok],
    scale = cf$scale[ok],
    shape = cf$shape,
    weights = cf$n[ok] / sum(cf$n[ok])
  )
  mixture_tail_quantile(mt, p)
}

# Transport: push each bin's conditional back through the copula, then combine
# the resulting marginal survival curves on the log scale.
est_transport <- function(cf, u, y, p) {
  ok <- is.finite(cf$threshold) &
    is.finite(cf$scale) &
    cf$tail_prob > 0 &
    cf$n > 0
  if (!any(ok) || !is.finite(cf$shape)) {
    return(rep(NA_real_, length(p)))
  }
  # Copula parameter from ranks only -- uses all n points, no tail involved.
  tau <- stats::cor(u, y, method = "kendall")
  rho_hat <- sin(pi * tau / 2)
  rho_hat <- max(min(rho_hat, 0.995), -0.995)

  idx <- which(ok)
  marginal_surv <- function(yy) {
    # Each bin gives its own estimate of the marginal survival; they agree in
    # expectation, so average their logs (a geometric mean). Averaging the
    # survivals themselves would let the heaviest bin dominate -- the same
    # failure mode the shared shape exists to avoid.
    ls <- vapply(
      idx,
      function(b) {
        s_cond <- cf$tail_prob[b] *
          gpd_survival(yy - cf$threshold[b], cf$scale[b], cf$shape)
        s_cond <- min(max(s_cond, 1e-300), 1 - 1e-12)
        log(max(transport_survival(s_cond, cf$u[b], rho_hat), 1e-300))
      },
      numeric(1L)
    )
    exp(mean(ls))
  }

  # Invert by bisection; the survival is monotone in yy.
  vapply(
    p,
    function(pp) {
      lo <- max(cf$threshold[ok], na.rm = TRUE)
      hi <- lo + 1
      it <- 0L
      while (marginal_surv(hi) > pp && it < 200L) {
        hi <- lo + (hi - lo) * 2
        it <- it + 1L
      }
      if (marginal_surv(lo) < pp) {
        return(NA_real_)
      }
      for (i in seq_len(120L)) {
        mid <- (lo + hi) / 2
        if (marginal_surv(mid) > pp) lo <- mid else hi <- mid
        if (hi - lo < 1e-8 * max(hi, 1)) break
      }
      (lo + hi) / 2
    },
    numeric(1L)
  )
}

# Anchored transport: take only the SHAPE from the copula decomposition, and
# anchor the level on the marginal body where data is plentiful.
#
# The full transport pushes an entire conditional curve through h^{-1}, so every
# error in the conditional fit is amplified by 1/CEVI -- a factor of 5 at
# rho = 0.9. But the marginal's moderate quantiles are very well determined
# empirically; it is only the index that is hard. So estimate the index by the
# decomposition, xi_Y = xi_cond / CEVI, and then fit the single remaining scale
# parameter to the marginal exceedances with that index held fixed.
#
# This is direct POT with the shape supplied by the copula instead of estimated
# from the handful of largest marginal observations.
est_anchored <- function(cf, u, y, p) {
  if (!is.finite(cf$shape)) {
    return(rep(NA_real_, length(p)))
  }
  tau <- stats::cor(u, y, method = "kendall")
  rho_hat <- max(min(sin(pi * tau / 2), 0.995), -0.995)
  cevi <- 1 - rho_hat^2
  if (!is.finite(cevi) || cevi < 0.02) {
    return(rep(NA_real_, length(p)))
  }
  xi_y <- cf$shape / cevi

  k <- max(10L, floor(POT_FRAC * length(y)))
  thr <- stats::quantile(y, 1 - k / length(y), names = FALSE)
  e <- y[y > thr] - thr
  zeta <- length(e) / length(y)
  sc <- gpd_profile_scale(e, rep(1, length(e)), xi_y)$scale
  if (!is.finite(sc)) {
    return(rep(NA_real_, length(p)))
  }
  thr + gpd_quantile_upper(p / zeta, sc, xi_y)
}

# --- run --------------------------------------------------------------------
true_level <- function(true_S, p) {
  vapply(
    p,
    function(pp) {
      lo <- 0
      hi <- 1
      while (true_S(hi) > pp) hi <- hi * 2
      for (i in seq_len(200L)) {
        mid <- (lo + hi) / 2
        if (true_S(mid) > pp) lo <- mid else hi <- mid
      }
      (lo + hi) / 2
    },
    numeric(1L)
  )
}

run_dgp <- function(dgp, label, note) {
  truth <- true_level(dgp(10L)$true_S, REPORT_P)
  methods <- c(
    "direct_pot",
    "dl_mixture",
    "transport",
    "transport+smooth",
    "anchored+smooth"
  )
  err <- array(NA_real_, c(N_REP, length(REPORT_P), length(methods)))
  shape_free <- shape_smooth <- numeric(N_REP)

  for (r in seq_len(N_REP)) {
    d <- dgp(N_OBS)
    cf <- fit_conditionals(d$u, d$y, scale_model = "free")
    cf_s <- fit_conditionals(d$u, d$y, scale_model = "loglinear")
    shape_free[r] <- cf$shape
    shape_smooth[r] <- cf_s$shape
    err[r, , 1] <- est_direct_pot(d$u, d$y, REPORT_P) / truth
    err[r, , 2] <- est_dl_mixture(cf, REPORT_P) / truth
    err[r, , 3] <- est_transport(cf, d$u, d$y, REPORT_P) / truth
    err[r, , 4] <- est_transport(cf_s, d$u, d$y, REPORT_P) / truth
    err[r, , 5] <- est_anchored(cf_s, d$u, d$y, REPORT_P) / truth
    if (r %% 50L == 0L) flush_cat(sprintf("   %s %d/%d\n", label, r, N_REP))
  }

  cat(sprintf("\n%s\n%s\n", label, note))
  cat(sprintf(
    "   true return levels: %s\n",
    paste(sprintf("%.1f", truth), collapse = "  ")
  ))
  cat(sprintf(
    "   fitted conditional EVI: free scales %.3f, log-linear %.3f\n\n",
    mean(shape_free, na.rm = TRUE),
    mean(shape_smooth, na.rm = TRUE)
  ))
  cat(sprintf("%-18s%s\n", "estimator", paste(
    sprintf("%26s", paste0("p = ", format(REPORT_P, scientific = TRUE))),
    collapse = ""
  )))
  cat(sprintf("%-18s%s\n", "", paste(
    rep(sprintf("%12s%14s", "median", "RMSE(log)"), length(REPORT_P)),
    collapse = ""
  )))
  cat(strrep("-", 18 + 26 * length(REPORT_P)), "\n")
  for (m in seq_along(methods)) {
    cells <- vapply(
      seq_along(REPORT_P),
      function(j) {
        e <- err[, j, m]
        e <- e[is.finite(e)]
        if (!length(e)) {
          return(sprintf("%12s%14s", "--", "--"))
        }
        sprintf("%11.2fx%14.2f", stats::median(e), sqrt(mean(log(e)^2)))
      },
      character(1L)
    )
    cat(sprintf("%-18s%s\n", methods[m], paste(cells, collapse = "")))
  }
  cat("\n   median = ratio to the true return level (1.00x is right).\n")
  cat("   RMSE(log) = root mean square of log(estimate/truth), so it is a\n")
  cat("   relative error and is comparable across return periods.\n")
}

if (!SWEEP_MODE) {
run_dgp(
  dgp_copula,
  "DGP A: Gaussian copula, exactly GP marginal",
  sprintf(
    paste(
      "   CEVI = %.2f, so conditionals are genuinely lighter than the marginal.",
      "   NOTE: the marginal here is exactly generalized Pareto, which is the",
      "   most favourable case direct POT could be handed.",
      sep = "\n"
    ),
    1 - RHO^2
  )
)

run_dgp(
  dgp_nuanced,
  "DGP C: Gaussian copula, marginal that is NOT close to GP",
  paste(
    "   Same CEVI, but the marginal is a two-component Pareto mixture whose",
    "   local EVI falls from 0.52 to 0.26 across the range a sample reaches,",
    "   approaching its asymptotic 0.25 from above. Direct POT no longer has",
    "   a GP target to hit.",
    sep = "\n"
  )
)

run_dgp(
  dgp_scale,
  "DGP B: conditional GP, bounded covariate effect (CEVI = 1)",
  paste(
    "   Marginal and conditional share an EVI, so the copula explains none of",
    "   the tail heaviness and transport has nothing to gain. This checks",
    "   whether it does harm when its premise fails.",
    sep = "\n"
  )
)

cat("=====================================================================\n")
}

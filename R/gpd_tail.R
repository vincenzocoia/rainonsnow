# Dependency-free generalized-Pareto kernels.
#
# These are the hot path for tail estimation and for evaluating mixture tails
# over many cells x many peak hours, so they deliberately avoid checkmate,
# distionary and the tidyverse: they take and return plain numeric vectors.
# The probaverse-facing wrappers live in mixture_tail.R.
#
# Parameterization matches distionary::dst_gp(): for excess y >= 0,
#   S(y) = (1 + xi * y / sigma) ^ (-1 / xi),           xi != 0
#   S(y) = exp(-y / sigma),                            xi == 0
# with support y in [0, Inf) when xi >= 0 and [0, -sigma/xi] when xi < 0.

# Treat |xi| below this as the exponential limit. The Gumbel-limit expressions
# lose relative precision well before xi underflows, so switch early.
XI_EPS <- 1e-8

#' Generalized Pareto survival function
#'
#' Vectorised over `y`, and over `scale` / `shape` by the usual recycling rules.
#'
#' @param y Numeric vector of excesses (values below `0` give survival `1`).
#' @param scale Positive numeric scale parameter.
#' @param shape Numeric shape parameter `xi`.
#' @returns Numeric vector of survival probabilities.
#' @examples
#' gpd_survival(c(0, 1, 10), scale = 2, shape = 0.2)
#' @export
gpd_survival <- function(y, scale, shape) {
  y <- pmax(y, 0)
  n <- max(length(y), length(scale), length(shape))
  y <- rep_len(y, n)
  scale <- rep_len(scale, n)
  shape <- rep_len(shape, n)

  out <- numeric(n)
  expo <- abs(shape) < XI_EPS
  if (any(expo)) {
    out[expo] <- exp(-y[expo] / scale[expo])
  }
  if (any(!expo)) {
    i <- !expo
    z <- 1 + shape[i] * y[i] / scale[i]
    # Past the upper endpoint (xi < 0) the survival is exactly zero.
    z <- pmax(z, 0)
    out[i] <- z^(-1 / shape[i])
  }
  out
}

#' Generalized Pareto density
#'
#' @inheritParams gpd_survival
#' @returns Numeric vector of densities (`0` outside the support).
#' @examples
#' gpd_density(c(0, 1, 10), scale = 2, shape = 0.2)
#' @export
gpd_density <- function(y, scale, shape) {
  n <- max(length(y), length(scale), length(shape))
  y <- rep_len(y, n)
  scale <- rep_len(scale, n)
  shape <- rep_len(shape, n)

  out <- numeric(n)
  inside <- y >= 0
  expo <- abs(shape) < XI_EPS

  i <- inside & expo
  if (any(i)) {
    out[i] <- exp(-y[i] / scale[i]) / scale[i]
  }
  i <- inside & !expo
  if (any(i)) {
    z <- 1 + shape[i] * y[i] / scale[i]
    ok <- z > 0
    val <- numeric(sum(i))
    val[ok] <- z[ok]^(-1 / shape[i][ok] - 1) / scale[i][ok]
    out[i] <- val
  }
  out
}

#' Generalized Pareto quantile from an upper-tail probability
#'
#' Returns `y` with `gpd_survival(y, scale, shape) == p`.
#'
#' @param p Numeric vector of upper-tail probabilities in `(0, 1]`.
#' @inheritParams gpd_survival
#' @returns Numeric vector of excesses.
#' @examples
#' gpd_quantile_upper(c(0.5, 0.01), scale = 2, shape = 0.2)
#' @export
gpd_quantile_upper <- function(p, scale, shape) {
  n <- max(length(p), length(scale), length(shape))
  p <- rep_len(p, n)
  scale <- rep_len(scale, n)
  shape <- rep_len(shape, n)

  out <- numeric(n)
  expo <- abs(shape) < XI_EPS
  if (any(expo)) {
    out[expo] <- -scale[expo] * log(p[expo])
  }
  if (any(!expo)) {
    i <- !expo
    out[i] <- scale[i] * (p[i]^(-shape[i]) - 1) / shape[i]
  }
  out
}

# ---------------------------------------------------------------------------
# Weighted likelihood
# ---------------------------------------------------------------------------

#' Weighted generalized Pareto negative log-likelihood
#'
#' Weights need not sum to one: they are the probability mass that an empirical
#' predictive distribution puts on each exceedance, so `sum(w)` is that
#' distribution's total tail mass.
#'
#' @param y Numeric vector of positive excesses.
#' @param w Numeric vector of non-negative weights, same length as `y`.
#' @param scale Positive scale parameter.
#' @param shape Shape parameter `xi`.
#' @returns A single numeric value; `Inf` when a point falls outside the
#'   support implied by `(scale, shape)`.
#' @examples
#' y <- c(0.2, 0.8, 1.5)
#' gpd_nllh(y, w = rep(1 / 3, 3), scale = 1, shape = 0.1)
#' @export
gpd_nllh <- function(y, w, scale, shape) {
  if (!is.finite(scale) || scale <= 0) {
    return(Inf)
  }
  if (abs(shape) < XI_EPS) {
    return(sum(w) * log(scale) + sum(w * y) / scale)
  }
  z <- 1 + shape * y / scale
  if (any(z <= 0)) {
    return(Inf)
  }
  sum(w) * log(scale) + (1 + 1 / shape) * sum(w * log(z))
}

# For a fixed shape, profile out the scale of a single component.
# The likelihood is well behaved in log(sigma), so a bounded 1-D search is both
# faster and more reliable here than a general-purpose optimiser.
#' @noRd
gpd_profile_scale <- function(y, w, shape, interval = NULL) {
  ybar <- sum(w * y) / sum(w)
  if (!is.finite(ybar) || ybar <= 0) {
    return(list(scale = NA_real_, nllh = Inf))
  }
  # sigma is on the order of the weighted mean excess; a factor of ~150 either
  # side comfortably brackets the optimum for any shape we allow.
  lo <- log(ybar) - 5
  hi <- log(ybar) + 5
  if (shape < 0) {
    # Support constraint: sigma > -shape * max(y).
    lo <- max(lo, log(-shape * max(y)) + 1e-6)
    if (lo >= hi) {
      hi <- lo + 10
    }
  }
  if (is.null(interval)) {
    interval <- c(lo, hi)
  }
  opt <- stats::optimize(
    function(ls) gpd_nllh(y, w, exp(ls), shape),
    interval = interval
  )
  list(scale = exp(opt$minimum), nllh = opt$objective)
}

#' Fit a generalized Pareto distribution to weighted excesses
#'
#' Independent (per-component) weighted maximum likelihood. This is the
#' single-distribution building block; [fit_gpd_shared_shape()] is the pooled
#' version that fixes the mixture tail problem described in
#' `vignette-free note` at the top of `R/mixture_tail.R`.
#'
#' @param y Numeric vector of positive excesses.
#' @param w Numeric vector of weights (defaults to equal weights).
#' @param shape_bounds Length-2 numeric search range for `xi`.
#' @returns A list with `scale`, `shape`, `nllh` and `n_eff`, the Kish
#'   effective sample size `sum(w)^2 / sum(w^2)` of the exceedances.
#' @examples
#' set.seed(1)
#' y <- gpd_quantile_upper(stats::runif(200), scale = 1, shape = 0.2)
#' fit_gpd_weighted(y)
#' @export
fit_gpd_weighted <- function(y, w = NULL, shape_bounds = c(-0.45, 1)) {
  keep <- is.finite(y) & y > 0
  if (is.null(w)) {
    w <- rep(1 / max(length(y), 1L), length(y))
  }
  keep <- keep & is.finite(w) & w > 0
  y <- y[keep]
  w <- w[keep]
  if (length(y) < 2L) {
    return(list(
      scale = NA_real_,
      shape = NA_real_,
      nllh = Inf,
      n_eff = 0
    ))
  }
  obj <- function(xi) gpd_profile_scale(y, w, xi)$nllh
  opt <- stats::optimize(obj, interval = shape_bounds)
  sc <- gpd_profile_scale(y, w, opt$minimum)
  list(
    scale = sc$scale,
    shape = opt$minimum,
    nllh = sc$nllh,
    n_eff = sum(w)^2 / sum(w^2)
  )
}

#' Fit generalized Pareto tails with a shape parameter shared across components
#'
#' Given several sets of weighted excesses -- in this project, the upper tail of
#' each peak hour's predictive distribution within one grid cell -- fit one
#' common shape `xi` and a separate scale for each component, by maximising the
#' summed weighted log-likelihood.
#'
#' The shape is profiled: for each candidate `xi` every component's scale is
#' optimised in closed-ish form by a bounded 1-D search, so the outer search is
#' only over the single parameter `xi`. That makes the cost linear in the number
#' of components rather than exponential in the parameter count, and it avoids
#' the flat, badly conditioned surfaces a joint `optim()` over `(xi, sigma_1,
#' ..., sigma_K)` runs into.
#'
#' @param excess_list List of numeric vectors of positive excesses, one per
#'   component.
#' @param weight_list List of numeric weight vectors matching `excess_list`.
#'   `NULL` means equal weights within each component.
#' @param shape_bounds Length-2 numeric search range for `xi`.
#' @param tol Tolerance for the outer 1-D search over `xi`.
#' @param bias_correct Whether to subtract the parametric-bootstrap estimate of
#'   the shape's bias. Defaults to `FALSE`; see the section below before
#'   turning it on.
#' @param n_boot Number of bootstrap replicates used when `bias_correct` is
#'   `TRUE`.
#'
#' @section Incidental parameters:
#' Because each component contributes its own scale, the number of nuisance
#' parameters grows with the number of components while the information per
#' component stays fixed. This is the Neyman--Scott situation, and the pooled
#' MLE of `xi` is biased low by O(1 / m) in the per-component effective sample
#' size `m` -- in simulation, about `-0.09` at `m = 15` and `-0.02` at
#' `m = 60`. Since a shape that is too small understates flood risk, the bias
#' is corrected by parametric bootstrap: simulate from the fitted model, refit,
#' and subtract the mean displacement. In simulation this brings the bias to
#' roughly `+0.005` with no increase in variance.
#'
#' A Cox--Reid adjusted profile likelihood was also tried for this and behaved
#' badly (the adjustment term dominated and drove `xi` to the search bound), so
#' the bootstrap is used instead.
#'
#' @section Why the correction is off by default:
#' Correcting the shape in isolation makes the *return levels* worse, which is
#' what the pipeline actually reports. There is a second, independent error:
#' each hour's scale is fitted from a handful of points, and the mixture return
#' level is convex in those scales, so scale noise alone inflates it (Jensen) by
#' about 7% at T = 100 years and 10% at T = 1000 -- measured in
#' `scripts/experiments/tail-index-pooling.R` as the `oracle_shape` row, where
#' the true shape is supplied and the level is still too high.
#'
#' The uncorrected pooled shape is biased low by roughly the same amount in the
#' opposite direction, so it lands close to correct: over T = 10, 100 and 1000
#' years it gives 1.00x, 0.99x and 0.96x of the true return level, against
#' 1.07x, 1.16x and 1.24x once the shape is corrected. Both are far better than
#' fitting a shape per hour, which reaches 1.64x at T = 1000 and keeps climbing.
#'
#' Relying on two errors cancelling is not satisfying, and the honest fix is to
#' correct the scale-noise inflation too -- most directly by bootstrapping the
#' return level itself rather than the shape. Until that exists, `FALSE` is the
#' better-calibrated default for return levels, and `TRUE` is the better choice
#' if the shape itself is the quantity of interest.
#'
#' @returns A list with:
#'   * `shape` -- the shared `xi` (bias-corrected when requested).
#'   * `shape_raw` -- the uncorrected pooled MLE.
#'   * `scale` -- numeric vector of per-component scales.
#'   * `nllh` -- the pooled negative log-likelihood at the optimum.
#'   * `n_eff` -- per-component effective sample sizes.
#'   * `shape_se` -- bootstrap standard error of `xi`, or `NA` when `n_boot` is
#'     zero. Reported for diagnostics; nothing in the pipeline consumes it.
#'   * `converged` -- `TRUE` when the optimum is interior to `shape_bounds`.
#' @examples
#' set.seed(1)
#' ex <- lapply(1:5, function(i) {
#'   gpd_quantile_upper(stats::runif(50), scale = i, shape = 0.15)
#' })
#' fit <- fit_gpd_shared_shape(ex, bias_correct = FALSE)
#' fit$shape
#' @export
fit_gpd_shared_shape <- function(
  excess_list,
  weight_list = NULL,
  shape_bounds = c(-0.45, 1),
  tol = 1e-4,
  bias_correct = FALSE,
  n_boot = 25L
) {
  if (is.null(weight_list)) {
    weight_list <- lapply(excess_list, function(y) {
      rep(1 / max(length(y), 1L), length(y))
    })
  }
  if (length(excess_list) != length(weight_list)) {
    stop("`excess_list` and `weight_list` must have the same length.")
  }

  # Drop components that cannot inform a two-parameter tail.
  clean <- Map(
    function(y, w) {
      keep <- is.finite(y) & y > 0 & is.finite(w) & w > 0
      list(y = y[keep], w = w[keep])
    },
    excess_list,
    weight_list
  )
  usable <- vapply(clean, function(e) length(e$y) >= 2L, logical(1L))
  clean_ok <- clean[usable]

  if (!length(clean_ok)) {
    return(list(
      shape = NA_real_,
      shape_raw = NA_real_,
      scale = rep(NA_real_, length(excess_list)),
      nllh = Inf,
      n_eff = rep(0, length(excess_list)),
      shape_se = NA_real_,
      converged = FALSE
    ))
  }

  profile <- function(xi) {
    total <- 0
    for (e in clean_ok) {
      total <- total + gpd_profile_scale(e$y, e$w, xi)$nllh
      if (!is.finite(total)) {
        return(Inf)
      }
    }
    total
  }

  opt <- stats::optimize(profile, interval = shape_bounds, tol = tol)
  xi <- opt$minimum

  scale <- rep(NA_real_, length(excess_list))
  scale[usable] <- vapply(
    clean_ok,
    function(e) gpd_profile_scale(e$y, e$w, xi)$scale,
    numeric(1L)
  )

  n_eff <- vapply(
    clean,
    function(e) if (length(e$w)) sum(e$w)^2 / sum(e$w^2) else 0,
    numeric(1L)
  )

  f0 <- opt$objective

  interior <- xi > shape_bounds[1] + 10 * tol &&
    xi < shape_bounds[2] - 10 * tol

  # Bias and spread both come from the same parametric-bootstrap replicates.
  #
  # There is deliberately no profile-deviance standard error here. The
  # components are weighted predictive distributions whose weights sum to one,
  # not independent observations, so the deviance is on a "one unit per
  # component" scale rather than the true information scale and a deviance
  # interval comes out badly miscalibrated (SE ~ 0.74 on data worth ~18000
  # points). The bootstrap spread is calibrated to the resampling scheme that
  # actually generated the components.
  shape_raw <- xi
  shape_se <- NA_real_
  if (n_boot > 0L && is.finite(xi)) {
    star <- shared_shape_boot(
      clean_ok,
      xi,
      scale[usable],
      shape_bounds,
      tol,
      n_boot
    )
    if (length(star)) {
      shape_se <- stats::sd(star)
      if (isTRUE(bias_correct)) {
        xi <- xi - (mean(star) - xi)
      }
    }
  }

  list(
    shape = xi,
    shape_raw = shape_raw,
    scale = scale,
    nllh = f0,
    n_eff = n_eff,
    shape_se = shape_se,
    converged = interior
  )
}

# Parametric-bootstrap replicates of the pooled shape estimator.
#' @noRd
shared_shape_boot <- function(clean_ok, xi, scale_ok, shape_bounds, tol, n_boot) {
  scale_ok <- scale_ok[is.finite(scale_ok)]
  if (!length(scale_ok) || !is.finite(xi)) {
    return(numeric(0))
  }
  sizes <- vapply(clean_ok, function(e) length(e$y), integer(1L))
  wts <- lapply(clean_ok, function(e) e$w)
  # Resample under the fitted model, preserving each component's sample size
  # and weights, then see where the estimator lands.
  star <- vapply(
    seq_len(n_boot),
    function(b) {
      sim <- Map(
        function(m, sg) gpd_quantile_upper(stats::runif(m), sg, xi),
        sizes,
        scale_ok
      )
      fit_gpd_shared_shape(
        sim,
        wts,
        shape_bounds = shape_bounds,
        tol = tol,
        bias_correct = FALSE,
        n_boot = 0L
      )$shape
    },
    numeric(1L)
  )
  star[is.finite(star)]
}

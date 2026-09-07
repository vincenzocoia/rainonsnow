# Closed-form upper tail of an equal-weight mixture of grafted GP predictions.
#
# WHY THIS EXISTS
#
# A cell's marginal runoff distribution is the equal-weight mixture of its
# peak-hour predictive distributions, each of which is empirical below its graft
# point and generalized Pareto above it. Two problems follow from evaluating
# that mixture as a probaverse object.
#
# 1. Tail index. If component i has shape xi_i, then
#      S_i(x) ~ zeta_i * (sigma_i / (xi_i * x)) ^ (1 / xi_i),
#    so the mixture survival is regularly varying with index
#    min_i (1 / xi_i) -- that is, the mixture's shape is `max(xi_i)`. The single
#    heaviest component sets the tail of the whole cell, and because each xi_i
#    is estimated from a small effective sample, `max_i xi_hat_i` is the maximum
#    of a few hundred noisy estimates and is badly inflated. See
#    `scripts/experiments/tail-index-pooling.R` for the size of the effect.
#
#    With a shared shape xi the problem disappears: every component is
#    regularly varying with the same index, and the mixture is too, with the
#    scales combining through the power mean in `mixture_tail_gpd_approx()`.
#
# 2. Cost. Evaluating the mixture through probaverse builds K distribution
#    objects and dispatches K generic calls per evaluation point. Above
#    `max(graft_of)` the mixture survival is an elementary function of the
#    stored `(graft_of, tail_prob, gp_scale, gp_shape)` quadruple, so it can be
#    evaluated as one matrix operation and inverted by bisection. That is the
#    difference between a return-level curve taking minutes and taking
#    milliseconds, and it needs no distribution objects on disk at all.

#' Build a closed-form mixture tail
#'
#' Describes the upper tail of `sum_i weights_i * F_i` where each `F_i` places
#' mass `tail_prob[i]` above `threshold[i]`, distributed as a generalized Pareto
#' with the given scale and shape.
#'
#' @param threshold Numeric vector of graft points `u_i`.
#' @param tail_prob Numeric vector of exceedance probabilities
#'   `P(X_i > u_i)`.
#' @param scale Numeric vector of GP scales.
#' @param shape Numeric vector of GP shapes, or a single value to share across
#'   all components (the recommended case -- see the note at the top of this
#'   file).
#' @param weights Mixture weights; default equal. Rescaled to sum to one.
#' @returns An object of class `mixture_tail`.
#' @examples
#' mt <- mixture_tail(
#'   threshold = c(10, 12, 11),
#'   tail_prob = c(0.5, 0.5, 0.5),
#'   scale = c(2, 3, 2.5),
#'   shape = 0.15
#' )
#' mixture_tail_survival(mt, c(20, 50, 100))
#' @export
mixture_tail <- function(
  threshold,
  tail_prob,
  scale,
  shape,
  weights = NULL
) {
  k <- length(threshold)
  if (k == 0L) {
    stop("`threshold` must have at least one element.")
  }
  tail_prob <- rep_len(tail_prob, k)
  scale <- rep_len(scale, k)
  shape <- rep_len(shape, k)
  if (is.null(weights)) {
    weights <- rep(1 / k, k)
  }
  weights <- rep_len(weights, k)

  keep <- is.finite(threshold) &
    is.finite(tail_prob) &
    is.finite(scale) &
    is.finite(shape) &
    is.finite(weights) &
    scale > 0 &
    tail_prob > 0 &
    weights > 0
  if (!any(keep)) {
    stop("No usable mixture components.")
  }

  weights <- weights[keep]
  threshold <- threshold[keep]
  tail_prob <- tail_prob[keep]
  scale <- scale[keep]
  shape <- shape[keep]
  weights <- weights / sum(weights)

  out <- list(
    threshold = threshold,
    tail_prob = tail_prob,
    scale = scale,
    shape = shape,
    weights = weights,
    n_dropped = sum(!keep)
  )

  # Shift every component up to the common base threshold once, at construction.
  #
  # For x >= u0, component i contributes
  #   zeta_i * (1 + xi*(x - u_i)/sigma_i) ^ (-1/xi)
  #     = zeta0_i * (1 + xi*(x - u0)/sigma0_i) ^ (-1/xi)
  # with sigma0_i = sigma_i + xi*(u0 - u_i) and
  #      zeta0_i = zeta_i * S_i(u0 - u_i).
  # Folding the mixture weights into zeta0 as well leaves the survival as a
  # single matrix power followed by one matrix-vector product, which is the
  # difference between this being usable inside an animation loop and not.
  u0 <- max(threshold)
  out$u0 <- u0
  out$shared_shape <- length(unique(shape)) == 1L
  if (out$shared_shape) {
    out$xi0 <- shape[[1L]]
    out$sigma0 <- scale + out$xi0 * (u0 - threshold)
    out$a0 <- weights * tail_prob * gpd_survival(u0 - threshold, scale, shape)
    ok0 <- is.finite(out$sigma0) & out$sigma0 > 0 & is.finite(out$a0)
    out$sigma0 <- out$sigma0[ok0]
    out$a0 <- out$a0[ok0]
  }

  structure(out, class = "mixture_tail")
}

#' @export
print.mixture_tail <- function(x, ...) {
  cat("<mixture_tail>\n")
  cat("  components: ", length(x$threshold), "\n", sep = "")
  sh <- unique(round(x$shape, 8))
  if (length(sh) == 1L) {
    cat("  shape:      ", format(sh), " (shared)\n", sep = "")
  } else {
    cat(
      "  shape:      ",
      format(min(x$shape)),
      " to ",
      format(max(x$shape)),
      " (NOT shared: mixture index is set by the max)\n",
      sep = ""
    )
  }
  cat("  valid above:", format(mixture_tail_lower(x)), "\n")
  invisible(x)
}

#' Lowest value at which a mixture tail is exact
#'
#' The closed form only describes the region above every component's graft
#' point. Below it the mixture still has empirical mass that is not represented
#' here, so evaluation is refused rather than silently extrapolated.
#'
#' @param mt A `mixture_tail`.
#' @returns A single numeric value.
#' @examples
#' mt <- mixture_tail(c(10, 12), c(0.5, 0.5), c(2, 3), 0.15)
#' mixture_tail_lower(mt)
#' @export
mixture_tail_lower <- function(mt) max(mt$threshold)

#' Largest exceedance probability a mixture tail can describe
#'
#' `mixture_tail_quantile()` cannot answer for probabilities above this, because
#' the corresponding quantile falls below [mixture_tail_lower()].
#'
#' @param mt A `mixture_tail`.
#' @returns A single numeric value.
#' @examples
#' mt <- mixture_tail(c(10, 12), c(0.5, 0.5), c(2, 3), 0.15)
#' mixture_tail_max_prob(mt)
#' @export
mixture_tail_max_prob <- function(mt) {
  mixture_tail_survival(mt, mixture_tail_lower(mt))
}

#' Survival function of a mixture tail
#'
#' @param mt A `mixture_tail`.
#' @param x Numeric vector of values.
#' @returns Numeric vector of exceedance probabilities. Values below
#'   [mixture_tail_lower()] return `NA` -- the mixture has empirical mass there
#'   that this object does not carry.
#' @examples
#' mt <- mixture_tail(c(10, 12), c(0.5, 0.5), c(2, 3), 0.15)
#' mixture_tail_survival(mt, c(20, 100))
#' @export
mixture_tail_survival <- function(mt, x) {
  lower <- mixture_tail_lower(mt)
  out <- rep(NA_real_, length(x))
  ok <- is.finite(x) & x >= lower
  if (!any(ok)) {
    return(out)
  }
  xo <- x[ok]

  if (isTRUE(mt$shared_shape)) {
    # Fast path: components already shifted to a common threshold, so this is
    # one elementwise power on an (points x components) matrix and one
    # matrix-vector product. No per-component dispatch, no recycling machinery.
    tt <- xo - mt$u0
    xi <- mt$xi0
    if (abs(xi) < XI_EPS) {
      m <- exp(-outer(tt, 1 / mt$sigma0))
    } else {
      m <- 1 + xi * outer(tt, 1 / mt$sigma0)
      m[m < 0] <- 0
      m <- m^(-1 / xi)
    }
    out[ok] <- as.vector(m %*% mt$a0)
    return(out)
  }

  # General path: shapes differ between components, so each needs its own
  # exponent. Kept for comparing against the unpooled status quo.
  ex <- outer(xo, mt$threshold, `-`)
  surv <- gpd_survival(
    as.vector(ex),
    rep(mt$scale, each = length(xo)),
    rep(mt$shape, each = length(xo))
  )
  dim(surv) <- dim(ex)
  out[ok] <- as.vector(surv %*% (mt$weights * mt$tail_prob))
  out
}

#' Density of a mixture tail
#'
#' @param mt A `mixture_tail`.
#' @param x Numeric vector of values.
#' @returns Numeric vector of densities; `NA` below [mixture_tail_lower()],
#'   where the mixture still has empirical mass this object does not carry.
#' @examples
#' mt <- mixture_tail(c(10, 12), c(0.5, 0.5), c(2, 3), 0.15)
#' mixture_tail_density(mt, c(20, 100))
#' @export
mixture_tail_density <- function(mt, x) {
  lower <- mixture_tail_lower(mt)
  out <- rep(NA_real_, length(x))
  ok <- is.finite(x) & x >= lower
  if (!any(ok)) {
    return(out)
  }
  xo <- x[ok]

  if (isTRUE(mt$shared_shape)) {
    dens <- gpd_density(
      rep(xo - mt$u0, times = length(mt$sigma0)),
      rep(mt$sigma0, each = length(xo)),
      mt$xi0
    )
    dim(dens) <- c(length(xo), length(mt$sigma0))
    out[ok] <- as.vector(dens %*% mt$a0)
    return(out)
  }

  ex <- outer(xo, mt$threshold, `-`)
  dens <- gpd_density(
    as.vector(ex),
    rep(mt$scale, each = length(xo)),
    rep(mt$shape, each = length(xo))
  )
  dim(dens) <- dim(ex)
  out[ok] <- as.vector(dens %*% (mt$weights * mt$tail_prob))
  out
}

#' Single-GPD approximation to a mixture tail
#'
#' When the shape is shared, the mixture is regularly varying with that same
#' shape and its asymptotic scale is a power mean of the component scales:
#' `sigma_eff = (sum_i w_i zeta_i sigma_i^(1/xi)) ^ xi / zeta_eff`. This is exact
#' in the limit and is used to seed the root-finding in
#' [mixture_tail_quantile()]; it is also a compact summary of a cell's tail.
#'
#' @param mt A `mixture_tail`.
#' @returns A list with `threshold`, `tail_prob`, `scale` and `shape`.
#' @examples
#' mt <- mixture_tail(c(10, 12), c(0.5, 0.5), c(2, 3), 0.15)
#' mixture_tail_gpd_approx(mt)
#' @export
mixture_tail_gpd_approx <- function(mt) {
  u <- mixture_tail_lower(mt)
  xi <- stats::weighted.mean(mt$shape, mt$weights)
  zeta_eff <- mixture_tail_survival(mt, u)
  # Match the mixture exactly at the base point, then take the asymptotic scale
  # implied by the power mean of component scales.
  if (abs(xi) < XI_EPS) {
    sc <- stats::weighted.mean(mt$scale, mt$weights * mt$tail_prob)
  } else {
    # Shift each component to the common threshold before combining.
    sig_u <- mt$scale + mt$shape * (u - mt$threshold)
    sig_u <- pmax(sig_u, .Machine$double.eps)
    zeta_u <- mt$tail_prob *
      gpd_survival(u - mt$threshold, mt$scale, mt$shape)
    a <- sum(mt$weights * zeta_u * sig_u^(1 / xi))
    sc <- (a / zeta_eff)^xi
  }
  list(threshold = u, tail_prob = zeta_eff, scale = sc, shape = xi)
}

#' Compress a mixture tail to fewer components
#'
#' With a shared shape the tail depends on the components only through the
#' weighted distribution of their scales, so many thousands of peak hours can be
#' summarised by a few dozen representative components with no practical loss.
#' Cost is linear in the number of components, so this is what makes a
#' return-period animation interactive.
#'
#' Components are grouped into equal-count bins of `sigma0`. Each bin keeps the
#' total weight of its members and takes the scale that reproduces their
#' combined asymptotic contribution, `(sum(a_i sigma_i^(1/xi)) / sum(a_i))^xi`,
#' so the far tail is preserved exactly and the near-tail error is bounded by
#' the within-bin spread of the scales.
#'
#' This is an approximation, so it is opt-in. Check it with
#' [mixture_tail_compression_error()] before relying on it.
#'
#' @param mt A `mixture_tail` with a shared shape.
#' @param n_bins Number of representative components to keep.
#' @returns A `mixture_tail` with at most `n_bins` components.
#' @examples
#' set.seed(1)
#' mt <- mixture_tail(
#'   threshold = stats::runif(2000, 8, 12),
#'   tail_prob = rep(0.5, 2000),
#'   scale = stats::runif(2000, 1, 6),
#'   shape = 0.15
#' )
#' small <- mixture_tail_compress(mt, 64)
#' length(small$threshold)
#' @export
mixture_tail_compress <- function(mt, n_bins = 64L) {
  if (!isTRUE(mt$shared_shape)) {
    stop(
      "Compression needs a shared shape. Components with different shapes do ",
      "not share a tail index, which is the problem a shared shape exists to ",
      "fix -- see the note at the top of R/mixture_tail.R."
    )
  }
  k <- length(mt$sigma0)
  if (k <= n_bins) {
    return(mt)
  }

  xi <- mt$xi0
  ord <- order(mt$sigma0)
  bin <- ceiling(seq_along(ord) * n_bins / k)
  sig <- mt$sigma0[ord]
  a <- mt$a0[ord]

  a_sum <- as.vector(tapply(a, bin, sum))
  if (abs(xi) < XI_EPS) {
    # Exponential limit: the matching statistic is the plain weighted mean.
    sig_rep <- as.vector(tapply(a * sig, bin, sum)) / a_sum
  } else {
    m1 <- as.vector(tapply(a * sig^(1 / xi), bin, sum))
    sig_rep <- (m1 / a_sum)^xi
  }

  # Rebuild at the common threshold: thresholds are all u0, so tail_prob
  # carries the weight and the weights themselves are uniform.
  mixture_tail(
    threshold = rep(mt$u0, length(a_sum)),
    tail_prob = a_sum * length(a_sum),
    scale = sig_rep,
    shape = xi,
    weights = rep(1 / length(a_sum), length(a_sum))
  )
}

#' Largest relative survival error introduced by compression
#'
#' @param mt The original `mixture_tail`.
#' @param compressed The result of [mixture_tail_compress()].
#' @param p Exceedance probabilities to check over.
#' @returns The maximum relative error in the return level over `p`.
#' @examples
#' set.seed(1)
#' mt <- mixture_tail(
#'   threshold = stats::runif(2000, 8, 12),
#'   tail_prob = rep(0.5, 2000),
#'   scale = stats::runif(2000, 1, 6),
#'   shape = 0.15
#' )
#' mixture_tail_compression_error(mt, mixture_tail_compress(mt, 64))
#' @export
mixture_tail_compression_error <- function(
  mt,
  compressed,
  p = 10^seq(-1, -5, length.out = 25)
) {
  p <- p[p <= min(mixture_tail_max_prob(mt), mixture_tail_max_prob(compressed))]
  if (!length(p)) {
    return(NA_real_)
  }
  full <- mixture_tail_quantile(mt, p)
  approx <- mixture_tail_quantile(compressed, p)
  max(abs(approx - full) / full, na.rm = TRUE)
}

#' Quantile of a mixture tail from an upper-tail probability
#'
#' Inverts [mixture_tail_survival()]. The survival function is strictly
#' decreasing, so the root is unique.
#'
#' The whole probability vector is solved at once by vectorised bisection: each
#' iteration evaluates the survival at every candidate simultaneously, so a
#' 500-point return curve costs a few dozen matrix operations rather than 500
#' separate one-dimensional root searches. The bracket is seeded from
#' [mixture_tail_gpd_approx()], which is exact in the limit and therefore
#' usually correct to within a factor of two straight away.
#'
#' @param mt A `mixture_tail`.
#' @param p Numeric vector of upper-tail probabilities.
#' @param tol Relative tolerance on the returned value.
#' @returns Numeric vector of quantiles; `NA` where `p` exceeds
#'   [mixture_tail_max_prob()] (the answer would lie below the graft region,
#'   where the mixture still has empirical mass this object does not carry).
#' @examples
#' mt <- mixture_tail(c(10, 12), c(0.5, 0.5), c(2, 3), 0.15)
#' mixture_tail_quantile(mt, c(0.1, 0.01, 0.001))
#' @export
mixture_tail_quantile <- function(mt, p, tol = 1e-10) {
  lower <- mixture_tail_lower(mt)
  pmax_ <- mixture_tail_max_prob(mt)
  approx <- mixture_tail_gpd_approx(mt)

  out <- rep(NA_real_, length(p))
  ok <- is.finite(p) & p > 0 & p <= pmax_
  if (!any(ok)) {
    return(out)
  }
  po <- p[ok]

  # Seed from the asymptotic GP, which is exact as p -> 0.
  seed <- gpd_quantile_upper(po / approx$tail_prob, approx$scale, approx$shape)
  seed[!is.finite(seed) | seed <= 0] <- 1

  # Work with the excess above `lower` so the bracket stays positive.
  lo <- rep(0, length(po))
  hi <- seed * 2

  # Push `hi` out until it brackets every root. The seed is close, so this
  # normally exits on the first check.
  for (i in seq_len(200L)) {
    need <- mixture_tail_survival(mt, lower + hi) > po
    if (!any(need, na.rm = TRUE)) {
      break
    }
    lo[need] <- hi[need]
    hi[need] <- hi[need] * 4
  }

  # Vectorised bisection: every bracket is halved on the same pass, so the cost
  # is a few dozen matrix operations for the whole curve rather than one root
  # search per point.
  for (i in seq_len(120L)) {
    mid <- (lo + hi) / 2
    high_side <- mixture_tail_survival(mt, lower + mid) > po
    high_side[is.na(high_side)] <- FALSE
    lo[high_side] <- mid[high_side]
    hi[!high_side] <- mid[!high_side]
    if (all(hi - lo <= tol * pmax(hi, 1))) {
      break
    }
  }

  out[ok] <- lower + (lo + hi) / 2
  out
}

#' Return levels of a mixture tail
#'
#' @param mt A `mixture_tail`.
#' @param at Numeric vector of return periods on the *event* axis, i.e. already
#'   multiplied by the number of POT events per year (as elsewhere in this
#'   repository -- see [enframe_at_events()]).
#' @returns Numeric vector of return levels.
#' @examples
#' mt <- mixture_tail(c(10, 12), c(0.5, 0.5), c(2, 3), 0.15)
#' # 100-year level at 3 events per year
#' mixture_tail_return_level(mt, at = 100 * 3)
#' @export
mixture_tail_return_level <- function(mt, at) {
  mixture_tail_quantile(mt, 1 / at)
}

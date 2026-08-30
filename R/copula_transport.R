# Transporting conditional distributions back to the marginal, through a copula.
#
# WHY
#
# The marginal of Y cannot be recovered by mixing the learned conditionals over
# the observed covariates. Every conditional shares the conditional tail index,
# so a finite mixture of them is regularly varying with THAT index -- lighter
# than the marginal's by the factor CEVI, the copula conditional extreme value
# index. The marginal only gets its full weight from covariate values beyond any
# sample, and machine learning cannot go there.
#
# The copula can. With U = F_X(X), V = F_Y(Y) and h(v|u) = dC(u,v)/du,
#
#     F_{Y|X=x}(y) = h( F_Y(y) | u ),      u = F_X(x),
#
# so inverting h in its first argument turns any ONE conditional into an
# estimate of the whole marginal:
#
#     F_Y(y) = h^{-1}( F_{Y|X=x}(y) | u ).
#
# COMBINING
#
# Every x gives its own estimate of the same marginal, so they should be
# combined as repeated measurements, not mixed. Mixing is the law of total
# probability and belongs to the other operation entirely; averaging estimates
# of one quantity is a different question, and there the pointwise MEDIAN is
# available. It is well defined as a distribution -- a pointwise median of
# monotone survival curves is itself monotone and lies in [0, 1] -- and, unlike
# the mean, it is not dragged around by whichever x happened to produce the
# heaviest tail. That is the same failure mode the shared shape exists to
# avoid, reappearing one level up.

#' Conditional copula function and its inverse
#'
#' `cop_h(v, u, family, par)` is `h(v | u) = dC(u, v) / du`, the conditional
#' distribution of `V` given `U = u`. `cop_h_inverse(w, u, family, par)` solves
#' `h(v | u) = w` for `v`.
#'
#' These work in probability space. In the far upper tail use [cop_h_surv()]
#' instead: composing these two there loses the answer to cancellation.
#'
#' @param v,w,u Numeric vectors in `(0, 1)`, recycled against each other.
#' @param family `"gaussian"`, `"clayton"` (lower tail dependence) or
#'   `"survival_clayton"` (Clayton rotated 180 degrees, so upper tail
#'   dependence -- usually the relevant one for flood peaks). The survival-space
#'   pair also accepts a function `f(s, u)` supplying the transform directly,
#'   which is how [td_transport()] feeds in a known copula that belongs to no
#'   family; `par` is then unused.
#' @param par Copula parameter: the correlation for `"gaussian"`, theta > 0 for
#'   the Clayton families.
#' @returns A numeric vector.
#' @examples
#' cop_h(0.9, 0.5, "gaussian", 0.7)
#' cop_h_inverse(cop_h(0.9, 0.5, "gaussian", 0.7), 0.5, "gaussian", 0.7)
#' @export
cop_h <- function(v, u, family = "gaussian", par = 0.5) {
  n <- max(length(v), length(u))
  v <- rep_len(v, n)
  u <- rep_len(u, n)
  v <- pmin(pmax(v, 1e-12), 1 - 1e-12)
  u <- pmin(pmax(u, 1e-12), 1 - 1e-12)
  switch(
    family,
    gaussian = stats::pnorm(
      (stats::qnorm(v) - par * stats::qnorm(u)) / sqrt(1 - par^2)
    ),
    clayton = clayton_h(v, u, par),
    # Rotating by 180 degrees swaps the roles of the corners: the survival
    # copula's upper tail is the Clayton lower tail.
    survival_clayton = 1 - clayton_h(1 - v, 1 - u, par),
    stop("Unknown copula family: ", family, call. = FALSE)
  )
}

#' @rdname cop_h
#' @export
cop_h_inverse <- function(w, u, family = "gaussian", par = 0.5) {
  n <- max(length(w), length(u))
  w <- rep_len(w, n)
  u <- rep_len(u, n)
  w <- pmin(pmax(w, 1e-15), 1 - 1e-15)
  u <- pmin(pmax(u, 1e-12), 1 - 1e-12)
  switch(
    family,
    gaussian = stats::pnorm(
      stats::qnorm(w) * sqrt(1 - par^2) + par * stats::qnorm(u)
    ),
    clayton = clayton_h_inv(w, u, par),
    survival_clayton = 1 - clayton_h_inv(1 - w, 1 - u, par),
    stop("Unknown copula family: ", family, call. = FALSE)
  )
}

#' @noRd
clayton_h <- function(v, u, theta) {
  b <- u^(-theta) + v^(-theta) - 1
  u^(-theta - 1) * b^(-1 - 1 / theta)
}

# P(V > v | U = u) for plain Clayton, in survival space.
# D = v^-theta - 1 = expm1(-theta * log1p(-s_v)), and
# h = (1 + D u^theta) ^ (-(1+theta)/theta), so 1 - h comes from expm1.
#' @noRd
clayton_surv <- function(s_v, u, theta) {
  d <- expm1(-theta * log1p(-s_v))
  -expm1(-((1 + theta) / theta) * log1p(d * u^theta))
}

#' @noRd
clayton_surv_inv <- function(s_w, u, theta) {
  d <- expm1(-theta / (1 + theta) * log1p(-s_w)) / u^theta
  -expm1(-log1p(d) / theta)
}

#' @noRd
clayton_h_inv <- function(w, u, theta) {
  # h = u^(-t-1) * b^(-(1+t)/t) with b = u^-t + v^-t - 1, so
  # b = (w * u^(t+1)) ^ (-t/(1+t)) and v follows.
  b <- (w * u^(theta + 1))^(-theta / (1 + theta))
  (b - u^(-theta) + 1)^(-1 / theta)
}

#' Conditional survival through the copula, and its inverse
#'
#' `cop_h_surv(s_v, u, ...)` is `P(V > v | U = u)` written in terms of
#' `s_v = P(V > v) = 1 - v`, and `cop_h_surv_inverse()` goes back the other way.
#'
#' These exist because [cop_h()] and [cop_h_inverse()] cannot be composed in the
#' far tail. Recovering a marginal as `1 - h^{-1}(1 - s)` forms `1 - s` for a
#' tiny `s`, which rounds to one and destroys the answer: at `s = 1e-10` the
#' round trip comes back 130 times too large, and `1 - s` has already underflowed
#' before the inverse is reached. Everything here stays in survival space, so no
#' small probability is ever formed by subtraction from one.
#'
#' All three families are exact in survival space. Plain Clayton needs a little
#' rearranging to get there: writing `D = v^-theta - 1`, the Clayton h-function
#' is exactly `(1 + D u^theta) ^ (-(1+theta)/theta)`, with the `u^-theta` having
#' cancelled algebraically rather than numerically. `D` and the final complement
#' then both come from `expm1`/`log1p` without ever forming a small number by
#' subtraction.
#'
#' @param s_v,s_w Numeric vectors of exceedance probabilities.
#' @param u Numeric vector of `F_X(x)` values.
#' @inheritParams cop_h
#' @returns A numeric vector of exceedance probabilities.
#' @examples
#' cop_h_surv_inverse(cop_h_surv(1e-8, 0.4, "gaussian", 0.7), 0.4, "gaussian", 0.7)
#' @export
cop_h_surv <- function(s_v, u, family = "gaussian", par = 0.5) {
  s_v <- pmin(pmax(s_v, 1e-300), 1 - 1e-15)
  u <- pmin(pmax(u, 1e-12), 1 - 1e-12)
  if (is.function(family)) {
    return(family(s_v, u))
  }
  switch(
    family,
    gaussian = stats::pnorm(
      (stats::qnorm(s_v, lower.tail = FALSE) - par * stats::qnorm(u)) /
        sqrt(1 - par^2),
      lower.tail = FALSE
    ),
    # P(V > v | U = u) for the survival copula is the Clayton h-function
    # evaluated at the reflected arguments, with no subtraction anywhere.
    survival_clayton = clayton_h(s_v, 1 - u, par),
    clayton = clayton_surv(s_v, u, par),
    stop("Unknown copula family: ", family, call. = FALSE)
  )
}

#' @rdname cop_h_surv
#' @export
cop_h_surv_inverse <- function(s_w, u, family = "gaussian", par = 0.5) {
  s_w <- pmin(pmax(s_w, 1e-300), 1 - 1e-15)
  u <- pmin(pmax(u, 1e-12), 1 - 1e-12)
  if (is.function(family)) {
    return(family(s_w, u))
  }
  switch(
    family,
    gaussian = stats::pnorm(
      stats::qnorm(s_w, lower.tail = FALSE) * sqrt(1 - par^2) +
        par * stats::qnorm(u),
      lower.tail = FALSE
    ),
    survival_clayton = clayton_h_inv(s_w, 1 - u, par),
    clayton = clayton_surv_inv(s_w, u, par),
    stop("Unknown copula family: ", family, call. = FALSE)
  )
}

#' Transport a conditional survival probability to the marginal
#'
#' Given `P(Y > y | X = x)` and the rank `u = F_X(x)`, returns `P(Y > y)`. This
#' is [cop_h_surv_inverse()] under the name that says what it is for.
#'
#' @param s_cond Numeric vector of conditional exceedance probabilities.
#' @param u Numeric vector of `F_X(x)` values.
#' @inheritParams cop_h
#' @returns A numeric vector of marginal exceedance probabilities.
#' @examples
#' transport_survival(0.01, u = 0.5, family = "gaussian", par = 0.7)
#' @export
transport_survival <- function(s_cond, u, family = "gaussian", par = 0.5) {
  cop_h_surv_inverse(s_cond, u, family, par)
}

#' An ensemble of marginal estimates, one per covariate value
#'
#' Each covariate value transports its own conditional tail back through the
#' copula, giving its own estimate of the whole marginal of `Y`. This holds all
#' of them and combines them.
#'
#' @param u Numeric vector of `F_X(x)` for each conditional.
#' @param threshold,tail_prob,scale Numeric vectors, one entry per conditional:
#'   the generalized Pareto threshold, the probability above it, and the scale.
#' @param shape Tail shape, shared across conditionals (length 1) or one each.
#' @param family,par Copula used for the transport; see [cop_h()].
#' @param weights Optional relative weights, used only by
#'   `combine = "weighted_mean"`. Defaults to equal.
#' @returns An object of class `marginal_ensemble`.
#' @examples
#' me <- marginal_ensemble(
#'   u = c(0.3, 0.5, 0.7),
#'   threshold = c(9, 10, 11),
#'   tail_prob = rep(0.3, 3),
#'   scale = c(2, 2.5, 3),
#'   shape = 0.1,
#'   family = "gaussian",
#'   par = 0.7
#' )
#' marginal_ensemble_survival(me, 50)
#' @export
marginal_ensemble <- function(
  u,
  threshold,
  tail_prob,
  scale,
  shape,
  family = "gaussian",
  par = 0.5,
  weights = NULL
) {
  k <- length(u)
  threshold <- rep_len(threshold, k)
  tail_prob <- rep_len(tail_prob, k)
  scale <- rep_len(scale, k)
  shape <- rep_len(shape, k)
  if (is.null(weights)) {
    weights <- rep(1, k)
  }
  weights <- rep_len(weights, k)

  keep <- is.finite(u) &
    u > 0 &
    u < 1 &
    is.finite(threshold) &
    is.finite(tail_prob) &
    tail_prob > 0 &
    is.finite(scale) &
    scale > 0 &
    is.finite(shape) &
    is.finite(weights) &
    weights > 0
  if (!any(keep)) {
    stop("No usable conditionals.", call. = FALSE)
  }

  structure(
    list(
      u = u[keep],
      threshold = threshold[keep],
      tail_prob = tail_prob[keep],
      scale = scale[keep],
      shape = shape[keep],
      weights = weights[keep] / sum(weights[keep]),
      family = family,
      par = par,
      n_dropped = sum(!keep)
    ),
    class = "marginal_ensemble"
  )
}

#' @export
print.marginal_ensemble <- function(x, ...) {
  cat("<marginal_ensemble>\n")
  cat("  conditionals: ", length(x$u), "\n", sep = "")
  cat(
    "  copula:       ",
    if (is.function(x$family)) {
      "supplied as a function"
    } else {
      paste0(x$family, " (par = ", format(x$par), ")")
    },
    "\n",
    sep = ""
  )
  sh <- unique(round(x$shape, 8))
  cat(
    "  conditional shape: ",
    if (length(sh) == 1L) format(sh) else
      paste0(format(min(x$shape)), " to ", format(max(x$shape))),
    "\n",
    sep = ""
  )
  cat("  valid above:  ", format(marginal_ensemble_lower(x)), "\n", sep = "")
  invisible(x)
}

#' @rdname marginal_ensemble
#' @param me A `marginal_ensemble`.
#' @export
marginal_ensemble_lower <- function(me) max(me$threshold)

#' Marginal survival estimates from an ensemble
#'
#' `marginal_ensemble_paths()` returns one estimate per conditional -- a matrix
#' of `length(y)` rows by one column per covariate value. That spread is a free
#' diagnostic: every column estimates the same marginal, so disagreement between
#' them measures how far the conditional model and the copula are from being
#' mutually consistent.
#'
#' @param me A `marginal_ensemble`.
#' @param y Numeric vector of values.
#' @returns A numeric matrix.
#' @examples
#' me <- marginal_ensemble(c(0.3, 0.7), c(9, 11), rep(0.3, 2), c(2, 3), 0.1)
#' dim(marginal_ensemble_paths(me, c(20, 50)))
#' @export
marginal_ensemble_paths <- function(me, y) {
  k <- length(me$u)
  out <- matrix(NA_real_, length(y), k)
  for (j in seq_len(k)) {
    s_cond <- me$tail_prob[j] *
      gpd_survival(y - me$threshold[j], me$scale[j], me$shape[j])
    out[, j] <- transport_survival(s_cond, me$u[j], me$family, me$par)
  }
  out
}

#' Combined marginal survival from an ensemble
#'
#' @inheritParams marginal_ensemble_paths
#' @param combine How to combine the per-covariate estimates.
#'   `"median"` (default) takes the pointwise median across covariate values.
#'   `"mean"` and `"weighted_mean"` take the pointwise arithmetic mean.
#'
#' @section Why the median:
#' These are repeated estimates of one quantity, not components of a mixture, so
#' combining them is an averaging problem rather than an application of the law
#' of total probability. The arithmetic mean is dragged upward by whichever
#' covariate value produced the heaviest tail -- the same domination that makes
#' a mixture of conditionals inherit `max(xi_i)`. The geometric mean fails the
#' other way, since `log S` runs to `-Inf` for a light tail; in simulation it
#' returned 0.65 of the true 1e-4 return level. The pointwise median is immune
#' to both, and a pointwise median of monotone survival curves is still a
#' monotone survival curve.
#'
#' @returns A numeric vector of exceedance probabilities.
#' @examples
#' me <- marginal_ensemble(c(0.3, 0.5, 0.7), c(9, 10, 11), rep(0.3, 3),
#'                         c(2, 2.5, 3), 0.1)
#' marginal_ensemble_survival(me, c(20, 50))
#' @export
marginal_ensemble_survival <- function(
  me,
  y,
  combine = c("median", "mean", "weighted_mean")
) {
  combine <- match.arg(combine)
  paths <- marginal_ensemble_paths(me, y)
  switch(
    combine,
    median = apply(paths, 1, stats::median, na.rm = TRUE),
    mean = rowMeans(paths, na.rm = TRUE),
    weighted_mean = as.vector(paths %*% me$weights)
  )
}

#' Combined marginal quantile from an ensemble
#'
#' Inverts [marginal_ensemble_survival()] by bisection over the whole
#' probability vector at once.
#'
#' @inheritParams marginal_ensemble_survival
#' @param p Numeric vector of upper-tail probabilities.
#' @param tol Relative tolerance on the returned value.
#' @returns A numeric vector of quantiles, `NA` where `p` is above what the
#'   ensemble can describe.
#' @examples
#' me <- marginal_ensemble(c(0.3, 0.5, 0.7), c(9, 10, 11), rep(0.3, 3),
#'                         c(2, 2.5, 3), 0.1)
#' marginal_ensemble_quantile(me, c(0.01, 1e-4))
#' @export
marginal_ensemble_quantile <- function(
  me,
  p,
  combine = c("median", "mean", "weighted_mean"),
  tol = 1e-10
) {
  combine <- match.arg(combine)
  lower <- marginal_ensemble_lower(me)
  s_at <- function(z) marginal_ensemble_survival(me, z, combine)

  out <- rep(NA_real_, length(p))
  ok <- is.finite(p) & p > 0 & p < s_at(lower)
  if (!any(ok)) {
    return(out)
  }
  po <- p[ok]

  lo <- rep(0, length(po))
  hi <- rep(max(1, lower), length(po))
  for (i in seq_len(200L)) {
    need <- s_at(lower + hi) > po
    if (!any(need, na.rm = TRUE)) {
      break
    }
    lo[need] <- hi[need]
    hi[need] <- hi[need] * 4
  }
  for (i in seq_len(200L)) {
    mid <- (lo + hi) / 2
    high <- s_at(lower + mid) > po
    high[is.na(high)] <- FALSE
    lo[high] <- mid[high]
    hi[!high] <- mid[!high]
    if (all(hi - lo <= tol * pmax(hi, 1))) {
      break
    }
  }
  out[ok] <- lower + (lo + hi) / 2
  out
}

#' Build a marginal ensemble from learned conditional tails
#'
#' Takes the upper tails of a set of learned conditional distributions -- the
#' output shape of [dl_tail_pieces()], one per covariate value -- and returns
#' the [marginal_ensemble()] they imply.
#'
#' By default one shape is fitted across all of them ([fit_gpd_shared_shape()]).
#' That is the same choice the pipeline makes within a cell, and for the same
#' reason: the conditionals are assumed to share a tail index and differ in
#' scale, so estimating a shape per covariate value only adds noise, which the
#' combination step then has to survive. Setting `shared_shape = FALSE` fits
#' each independently, which is what makes the spread across covariate values an
#' honest diagnostic rather than a flat line by construction.
#'
#' @param pieces List of tail pieces, each with `threshold`, `excess`,
#'   `weight` and `tail_prob`, as [dl_tail_pieces()] returns. `NULL` entries are
#'   dropped.
#' @param u Numeric vector of `F_X(x)`, one per entry of `pieces`.
#' @param shape Optional tail shape to impose. When given, only the scales are
#'   fitted.
#' @param shared_shape Fit one shape across all conditionals (the default) or
#'   one each.
#' @param weights Optional relative weights; defaults to each conditional's
#'   effective sample size.
#' @inheritParams marginal_ensemble
#' @returns A `marginal_ensemble`.
#' @examples
#' pieces <- lapply(1:3, function(i) {
#'   list(threshold = 10 * i, excess = stats::rexp(50, 1 / i),
#'        weight = rep(1 / 50, 50), tail_prob = 0.3, n_eff = 50)
#' })
#' transport_ensemble(pieces, u = c(0.3, 0.5, 0.7), par = 0.7)
#' @export
transport_ensemble <- function(
  pieces,
  u,
  family = "gaussian",
  par = 0.5,
  shape = NULL,
  shared_shape = TRUE,
  weights = NULL
) {
  if (length(pieces) != length(u)) {
    stop("`pieces` and `u` must have the same length.", call. = FALSE)
  }
  keep <- !vapply(pieces, is.null, logical(1L))
  pieces <- pieces[keep]
  u <- u[keep]
  if (!length(pieces)) {
    stop("No usable conditionals.", call. = FALSE)
  }
  if (is.null(weights)) {
    weights <- vapply(
      pieces,
      function(p) {
        if (is.null(p$n_eff)) sum(p$weight)^2 / sum(p$weight^2) else p$n_eff
      },
      numeric(1L)
    )
  } else {
    weights <- rep_len(weights, length(keep))[keep]
  }

  excess <- lapply(pieces, `[[`, "excess")
  weight <- lapply(pieces, `[[`, "weight")

  if (!is.null(shape)) {
    scale <- vapply(
      seq_along(pieces),
      function(i) gpd_profile_scale(excess[[i]], weight[[i]], shape)$scale,
      numeric(1L)
    )
    shape <- rep(shape, length(pieces))
  } else if (shared_shape) {
    fit <- fit_gpd_shared_shape(excess, weight, n_boot = 0L)
    scale <- fit$scale
    shape <- rep(fit$shape, length(pieces))
  } else {
    fits <- lapply(
      seq_along(pieces),
      function(i) fit_gpd_weighted(excess[[i]], weight[[i]])
    )
    scale <- vapply(fits, `[[`, numeric(1L), "scale")
    shape <- vapply(fits, `[[`, numeric(1L), "shape")
  }

  marginal_ensemble(
    u = u,
    threshold = vapply(pieces, `[[`, numeric(1L), "threshold"),
    tail_prob = vapply(pieces, `[[`, numeric(1L), "tail_prob"),
    scale = scale,
    shape = shape,
    family = family,
    par = par,
    weights = weights
  )
}

#' Copula density
#'
#' Needed to fit a copula by likelihood on part of the unit square rather than
#' from a body statistic like Kendall's tau. See [fit_copula_upper()].
#'
#' @param u,v Numeric vectors in `(0, 1)`, recycled against each other.
#' @inheritParams cop_h
#' @returns A numeric vector.
#' @examples
#' cop_density(0.9, 0.95, "gaussian", 0.7)
#' @export
cop_density <- function(u, v, family = "gaussian", par = 0.5) {
  u <- pmin(pmax(u, 1e-12), 1 - 1e-12)
  v <- pmin(pmax(v, 1e-12), 1 - 1e-12)
  switch(
    family,
    gaussian = {
      a <- stats::qnorm(u)
      b <- stats::qnorm(v)
      r <- par
      exp(
        -(r^2 * (a^2 + b^2) - 2 * r * a * b) / (2 * (1 - r^2)) -
          0.5 * base::log(1 - r^2)
      )
    },
    clayton = clayton_density(u, v, par),
    survival_clayton = clayton_density(1 - u, 1 - v, par),
    stop("Unknown copula family: ", family, call. = FALSE)
  )
}

#' @noRd
clayton_density <- function(u, v, theta) {
  (1 + theta) *
    (u * v)^(-1 - theta) *
    (u^-theta + v^-theta - 1)^(-1 / theta - 2)
}

#' Fit a copula parameter on the upper region only
#'
#' Kendall's tau is a statistic of the whole unit square, and the transport
#' depends on the copula near the top. This fits the parameter by maximum
#' likelihood restricted to observations with `v` above a threshold, correctly
#' normalised for that region: the density of `V` given `U = u` and `V > v0` is
#' `c(u, v) / P(V > v0 | U = u)`.
#'
#' The restriction is on `v` alone, across all of `u`. That is deliberate. The
#' upper-right corner is not the only place the copula shapes the conditional
#' tail index -- for families with a non-constant conditional extreme value
#' index the relevant behaviour is spread along the upper edge (Coia, Joe and
#' Nolde 2024) -- so conditioning on a joint corner would throw away the part of
#' the sample that carries it.
#'
#' @param u,v Pseudo-observations in `(0, 1)`.
#' @param family Copula family; see [cop_h()].
#' @param v_threshold Keep observations with `v` above this.
#' @param bounds Search interval for the parameter.
#' @returns A single number, `NA` if too few observations remain.
#' @examples
#' set.seed(1)
#' u <- stats::runif(300)
#' v <- cop_h_inverse(stats::runif(300), u, "gaussian", 0.6)
#' fit_copula_upper(u, v, "gaussian")
#' @export
fit_copula_upper <- function(
  u,
  v,
  family = "gaussian",
  v_threshold = 0.5,
  bounds = NULL
) {
  keep <- is.finite(u) & is.finite(v) & v > v_threshold
  if (sum(keep) < 8L) {
    return(NA_real_)
  }
  uu <- u[keep]
  vv <- v[keep]
  if (is.null(bounds)) {
    bounds <- if (identical(family, "gaussian")) c(-0.98, 0.98) else c(0.02, 40)
  }
  nllh <- function(par) {
    d <- cop_density(uu, vv, family, par)
    norm <- cop_h_surv(1 - v_threshold, uu, family, par)
    if (any(!is.finite(d)) || any(d <= 0) || any(norm <= 0)) {
      return(1e10)
    }
    -sum(base::log(d) - base::log(norm))
  }
  o <- try(stats::optimize(nllh, interval = bounds), silent = TRUE)
  if (inherits(o, "try-error")) {
    return(NA_real_)
  }
  o$minimum
}

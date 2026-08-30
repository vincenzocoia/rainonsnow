# A rain-on-snow shaped process: two drivers, one of them bounded.
#
# THE PHYSICAL STORY, AND WHY IT MATTERS FOR THE TAIL
#
# An annual peak flow is the larger of two things that can cause it. A
# snowmelt-driven peak is limited by how much snow there is to melt, so its
# distribution has a finite upper endpoint. A rain-driven peak has no such
# ceiling and inherits the heavy tail of the rainfall that caused it. Writing X
# for the rainfall amount, which is observed,
#
#     A ~ GP(scale_a, shape_a),  shape_a < 0     snowmelt-driven, BOUNDED
#     B | X = x ~ GP(b_coef * x, shape_b)        rain-driven, scale set by rain
#     Y = max(A, B)                              the annual peak
#
# with A independent of X. Everything else is derived: the conditional survival
# is exact and elementary,
#
#     P(Y > y | X = x) = 1 - F_A(y) F_{B|X}(y | x),
#
# and the marginal is that integrated over X, by quadrature.
#
# WHAT THIS CONSTRUCTION IS FOR
#
# The marginal of Y has a kink. Below A's endpoint it is dominated by the
# bounded driver and looks light; above it, only the rain-driven component
# survives and it is heavy. A generalized Pareto fitted to Y over a range that
# straddles the kink reads a shape somewhere in between and then extrapolates
# it, which understates the long return periods. At the settings used in
# `scripts/experiments/hydrological-decomposition.R` the asymptotic bias of a
# top-25% fit is +16% at the 20-year level and -24% at the 200-year level: it is
# not a variance problem and more data does not remove it.
#
# The covariate does something a generic copula cannot. It selects the REGIME.
# At low rainfall the conditional is snowmelt-dominated and its best generalized
# Pareto shape is negative; at high rainfall it is rain-dominated and the shape
# is positive. Each conditional is a single clean process where the marginal is
# a mixture of two, which is the decomposition argument in its strongest form --
# and note that it needs the conditional SHAPE to move with the covariate, not
# just the scale.
#
# Reference for the tail-index bookkeeping: Coia, Joe and Nolde (2024),
# Copula-based conditional tail indices, J. Multivariate Analysis 201, 105268.

#' A two-driver process with one bounded driver
#'
#' @param scale_a,shape_a Generalized Pareto parameters of the snowmelt-driven
#'   peak. A negative `shape_a` gives it the finite upper endpoint
#'   `scale_a / |shape_a|`, which is what creates the kink.
#' @param scale_x,shape_x Generalized Pareto parameters of the rainfall `X`.
#' @param b_coef,shape_b The rain-driven peak given rainfall `x` is generalized
#'   Pareto with scale `b_coef * x` and shape `shape_b`.
#' @returns An object of class `two_process_dgp`.
#' @examples
#' d <- two_process_dgp()
#' tp_marginal_quantile(d, 1 / c(20, 100))
#' @export
two_process_dgp <- function(
  scale_a = 1,
  shape_a = -0.25,
  scale_x = 1,
  shape_x = 0.25,
  b_coef = 0.5,
  shape_b = 0.10
) {
  if (shape_a >= 0) {
    stop(
      "`shape_a` must be negative: the bounded driver is the point.",
      call. = FALSE
    )
  }
  if (scale_a <= 0 || scale_x <= 0 || b_coef <= 0) {
    stop("Scales must be positive.", call. = FALSE)
  }
  dgp <- structure(
    list(
      scale_a = scale_a,
      shape_a = shape_a,
      scale_x = scale_x,
      shape_x = shape_x,
      b_coef = b_coef,
      shape_b = shape_b,
      endpoint_a = scale_a / abs(shape_a)
    ),
    class = "two_process_dgp"
  )
  # Cache log S_Y on a log-y grid: the quadrature is exact but too slow to sit
  # inside a bisection. `exact = TRUE` still runs it, and the tests check they
  # agree.
  ly <- seq(log(1e-4), log(1e9), length.out = 3000L)
  dgp$cache <- list(
    log_y = ly,
    log_s = tp_marginal_survival(dgp, exp(ly), log = TRUE, exact = TRUE)
  )
  dgp
}

#' @export
print.two_process_dgp <- function(x, ...) {
  cat("<two_process_dgp>\n")
  cat(sprintf(
    "  A (snowmelt-driven) ~ GP(%g, %g), bounded above by %.3f\n",
    x$scale_a,
    x$shape_a,
    x$endpoint_a
  ))
  cat(sprintf("  X (rainfall)        ~ GP(%g, %g)\n", x$scale_x, x$shape_x))
  cat(sprintf(
    "  B | X = x           ~ GP(%g * x, %g)\n",
    x$b_coef,
    x$shape_b
  ))
  cat(sprintf(
    "  Y = max(A, B); marginal tail index %.3f\n",
    max(x$shape_x, x$shape_b)
  ))
  invisible(x)
}

#' The rainfall distribution
#'
#' @param dgp A `two_process_dgp`.
#' @param p Numeric vector of probabilities.
#' @param x Numeric vector of rainfall values.
#' @returns A numeric vector.
#' @examples
#' tp_x_quantile(two_process_dgp(), c(0.1, 0.9))
#' @export
tp_x_quantile <- function(dgp, p) {
  gpd_quantile_upper(1 - p, dgp$scale_x, dgp$shape_x)
}

#' @rdname tp_x_quantile
#' @export
tp_x_cdf <- function(dgp, x) 1 - gpd_survival(x, dgp$scale_x, dgp$shape_x)

#' Conditional survival of the annual peak given rainfall
#'
#' `P(Y > y | X = x) = 1 - F_A(y) F_{B|X}(y | x)`, exact.
#'
#' @inheritParams tp_x_quantile
#' @param y,x Numeric vectors, recycled against each other.
#' @returns A numeric vector.
#' @examples
#' tp_conditional_survival(two_process_dgp(), y = c(2, 5), x = 1)
#' @export
tp_conditional_survival <- function(dgp, y, x) {
  sa <- gpd_survival(y, dgp$scale_a, dgp$shape_a)
  sb <- gpd_survival(y, dgp$b_coef * x, dgp$shape_b)
  sa + sb - sa * sb
}

#' Marginal survival of the annual peak
#'
#' The conditional integrated over the rainfall distribution, by quadrature in
#' `log X`. Nothing about this marginal is chosen; it is what the two drivers
#' imply, and it belongs to no standard family.
#'
#' @inheritParams tp_conditional_survival
#' @param ngrid Number of quadrature points.
#' @param log Return `log S_Y(y)` rather than `S_Y(y)`.
#' @param exact Run the quadrature instead of interpolating the cached curve.
#' @returns A numeric vector.
#' @examples
#' tp_marginal_survival(two_process_dgp(), c(1, 5, 20))
#' @export
tp_marginal_survival <- function(
  dgp,
  y,
  ngrid = 8001L,
  log = FALSE,
  exact = FALSE
) {
  if (!exact && !is.null(dgp$cache)) {
    out <- stats::approx(
      dgp$cache$log_y,
      dgp$cache$log_s,
      xout = base::log(y),
      rule = 2
    )$y
    out[base::log(y) < dgp$cache$log_y[1L]] <- 0
    return(if (log) out else exp(out))
  }
  l <- seq(-14, 60, length.out = ngrid)
  dl <- l[2L] - l[1L]
  # Density of L = log X, which is f_X(e^l) e^l.
  lf <- base::log(gpd_density(exp(l), dgp$scale_x, dgp$shape_x)) + l
  ok <- is.finite(lf)
  l <- l[ok]
  lf <- lf[ok]
  w <- rep(dl, length(l))
  w[c(1L, length(w))] <- dl / 2
  fx <- exp(lf) * w
  total <- sum(fx)
  x <- exp(l)
  out <- vapply(
    y,
    function(yy) base::log(sum(fx * tp_conditional_survival(dgp, yy, x)) / total),
    numeric(1L)
  )
  if (log) out else exp(out)
}

#' Marginal quantile of the annual peak
#'
#' @inheritParams tp_marginal_survival
#' @param s Numeric vector of upper-tail probabilities. For annual maxima the
#'   `T`-year level is `s = 1 / T`.
#' @returns A numeric vector.
#' @examples
#' tp_marginal_quantile(two_process_dgp(), 1 / c(20, 100, 200))
#' @export
tp_marginal_quantile <- function(dgp, s) {
  exp(stats::approx(
    rev(dgp$cache$log_s),
    rev(dgp$cache$log_y),
    xout = base::log(s),
    rule = 2
  )$y)
}

#' Local tail index of the marginal
#'
#' @inheritParams tp_marginal_survival
#' @param eps Half-width of the central difference, in `log y`.
#' @returns A numeric vector.
#' @examples
#' td <- two_process_dgp()
#' tp_local_evi(td, tp_marginal_quantile(td, c(0.05, 0.005)))
#' @export
tp_local_evi <- function(dgp, y, eps = 0.02) {
  lo <- tp_marginal_survival(dgp, y * exp(-eps), log = TRUE)
  hi <- tp_marginal_survival(dgp, y * exp(eps), log = TRUE)
  2 * eps / (lo - hi)
}

#' Simulate from the process
#'
#' @inheritParams tp_x_quantile
#' @param n Number of years.
#' @returns A data frame with the rainfall `x`, both drivers `a` and `b`, the
#'   annual peak `y`, and `u = F_X(x)`. The drivers are returned because one
#'   experiment treats them as observed.
#' @examples
#' head(tp_simulate(two_process_dgp(), 5))
#' @export
tp_simulate <- function(dgp, n) {
  u <- stats::runif(n)
  x <- tp_x_quantile(dgp, u)
  a <- gpd_quantile_upper(stats::runif(n), dgp$scale_a, dgp$shape_a)
  b <- gpd_quantile_upper(stats::runif(n), dgp$b_coef * x, dgp$shape_b)
  data.frame(x = x, a = a, b = b, y = pmax(a, b), u = u)
}

#' The process's own copula, in survival space
#'
#' Takes `P(Y > y | X = x)` and `u = F_X(x)` and returns `P(Y > y)`, exactly,
#' by undoing the conditional and evaluating the derived marginal. Pass this to
#' [marginal_ensemble()] as `family` to transport under the true dependence.
#'
#' @inheritParams tp_conditional_survival
#' @param s_cond Numeric vector of conditional exceedance probabilities.
#' @param u Numeric vector of `F_X(x)`.
#' @returns A numeric vector of marginal exceedance probabilities.
#' @examples
#' d <- two_process_dgp()
#' tp_transport(d, 0.01, u = 0.9)
#' @export
tp_transport <- function(dgp, s_cond, u) {
  n <- max(length(s_cond), length(u))
  s_cond <- rep_len(s_cond, n)
  u <- rep_len(u, n)
  x <- tp_x_quantile(dgp, u)
  y <- vapply(
    seq_len(n),
    function(i) {
      f <- function(z) tp_conditional_survival(dgp, exp(z), x[i]) - s_cond[i]
      lo <- -30
      hi <- 40
      if (f(lo) < 0) {
        return(0)
      }
      if (f(hi) > 0) {
        return(Inf)
      }
      for (k in seq_len(200L)) {
        m <- (lo + hi) / 2
        if (f(m) > 0) lo <- m else hi <- m
        if (hi - lo < 1e-10) {
          break
        }
      }
      exp((lo + hi) / 2)
    },
    numeric(1L)
  )
  tp_marginal_survival(dgp, y)
}

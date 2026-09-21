#' Convert an ecdf to an empirical distribution
#'
#' @param fun An ecdf function.
#' @returns An empirical distribution from the distionary package.
#' @examples
#' fun <- ecdf(rnorm(100))
#' ecdf_to_dst(fun)
#' @export
ecdf_to_dst <- function(fun) {
  checkmate::assert_class(fun, "ecdf")
  k <- stats::knots(fun)
  fun_k <- fun(k)
  steps0 <- diff(fun_k)
  steps <- append(fun_k[1], steps0)
  distionary::dst_empirical(k, weights = steps)
}

#' Check if rows of a data frame have missing values
#'
#' @param df A data frame.
#' @returns A logical vector of length `nrow(df)` indicating which rows have
#'   missing values.
#' @examples
#' df <- data.frame(x = c(1, 2, NA, 4), y = c(NA, 2, 3, 4))
#' df_rows_have_missing(df)
#' @export
df_rows_have_missing <- function(df) {
  lgl_mat <- as.matrix(dplyr::mutate(
    df,
    dplyr::across(tidyselect::everything(), is.na)
  ))
  apply(lgl_mat, 1, any)
}

#' Return periods for reporting
#'
#' Canonical year-based grid used in summaries and maps. This is the default
#' for [enframe_at_events()] (`return_periods` argument): each value is multiplied
#' by `num_events_per_year` there to evaluate return levels on a POT event axis.
#'
#' @returns A numeric vector of commonly reported return periods (years).
#' @seealso [enframe_at_events()]
#' @examples
#' rp_reporting()
#' @export
rp_reporting <- function() c(2, 5, 10, 20, 50, 100, 200, 500, 1000)

#' Return periods for marginal frequency–magnitude curves
#'
#' Log-spaced grid from 2 to 500 years (500 points), used when precomputing or
#' plotting DL and naive marginal return curves in
#' `scripts/5-runoff_marginals.r` and `apps/5-runoff-marginals-explorer`.
#' Maps and summaries should still use [rp_reporting()].
#'
#' @returns A numeric vector of return periods in years.
#' @seealso [rp_reporting()], [enframe_at_events()]
#' @examples
#' rp_marginal_curve()
#' length(rp_marginal_curve())
#' @export
rp_marginal_curve <- function() {
  exp(seq(log(2), log(500), length.out = 500))
}

#' Return levels by numerical CDF inversion
#'
#' Drop-in replacement for `distionary::eval_return()` for distributions whose
#' quantile function is disabled in the installed distionary (non-continuous
#' mixtures "lacking structured support", e.g. the GP-tail marginal built by
#' [mix2()] on grafted distributions). A return level for period `r` is the
#' quantile at `p = 1 - 1/r`, found here by root-finding on the CDF (which
#' remains available) rather than the disabled quantile algorithm. Matches
#' `eval_return()` to numerical precision in the tail, where return levels lie.
#'
#' @param distribution A distionary distribution (`dst`) with a working CDF.
#' @param at Numeric vector of return periods (same convention as
#'   `distionary::eval_return()`; the event-frequency axis in
#'   `scripts/5-runoff_marginals.r`).
#' @param lower Lower bracket for the root search (default `0`, appropriate for
#'   non-negative runoff).
#' @param n_grid Number of grid points used to tabulate the CDF before
#'   inverting it by interpolation.
#' @returns A numeric vector of return levels, one per element of `at`.
#' @seealso [rp_marginal_curve()]
#' @export
eval_return_numeric <- function(distribution, at, lower = 0, n_grid = 4000L) {
  checkmate::assert_class(distribution, "dst")
  checkmate::assert_numeric(at)
  p <- 1 - 1 / at
  cdf_at <- function(x) as.numeric(distionary::eval_cdf(distribution, at = x))

  # Grow an upper bound until the CDF covers the largest requested probability,
  # then evaluate the (expensive) mixture CDF once on a dense grid and invert by
  # interpolation. One vectorised CDF call per distribution instead of thousands.
  max_p <- suppressWarnings(max(p[is.finite(p) & p < 1], na.rm = TRUE))
  hi <- max(1, abs(lower) + 1)
  iter <- 0L
  while (is.finite(max_p) && cdf_at(hi) < max_p && iter < 80L) {
    hi <- hi * 2
    iter <- iter + 1L
  }

  x_grid <- seq(lower, hi, length.out = n_grid)
  cdf_grid <- cdf_at(x_grid)
  # Keep a strictly increasing CDF so interpolation is well defined (flat regions
  # of a step/mixture CDF create ties).
  keep <- c(TRUE, diff(cdf_grid) > 0)
  x_grid <- x_grid[keep]
  cdf_grid <- cdf_grid[keep]

  out <- rep(NA_real_, length(p))
  ok <- is.finite(p) & p > 0 & p < 1
  out[is.finite(p) & p <= 0] <- lower
  out[is.finite(p) & p >= 1] <- Inf
  if (any(ok) && length(cdf_grid) >= 2) {
    out[ok] <- stats::approx(
      x = cdf_grid,
      y = x_grid,
      xout = p[ok],
      rule = 2
    )$y
  }
  out
}

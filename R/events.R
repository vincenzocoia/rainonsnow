#' T-year runoff event levels from a direct POT fit
#'
#' Fits a generalized Pareto distribution to one cell's peak excesses over its
#' POT threshold and converts return periods (years) to runoff levels using the
#' peak rate. This defines "the T-year event" independently of any
#' distributional learning model, so that every analysis queries the same
#' event and results can be compared across analyses.
#'
#' With peaks arriving at `rate` per year and excess distribution `G`, the
#' T-year level `z_T` solves `rate * (1 - G(z_T - u)) = 1 / T`.
#'
#' @param runoff Peak runoff values for one cell (all above `threshold`).
#' @param threshold The POT threshold `u` for the cell.
#' @param n_years Length of the record in years.
#' @param return_periods Return periods in years.
#' @returns A tibble with `return_period`, `return_level`, `rate`, `gp_scale`,
#'   `gp_shape`.
#' @export
pot_event_levels <- function(runoff, threshold, n_years, return_periods) {
  excess <- runoff - threshold
  excess <- excess[is.finite(excess) & excess >= 0]
  gp <- famish::fit_dst_gp(excess, method = "mle")
  rate <- length(excess) / n_years
  p_exceed <- 1 / (rate * return_periods)
  lvl <- rep(NA_real_, length(return_periods))
  ok <- p_exceed < 1
  lvl[ok] <- threshold + distionary::eval_quantile(gp, at = 1 - p_exceed[ok])
  pars <- distionary::parameters(gp)
  tibble::tibble(
    return_period = return_periods,
    return_level = lvl,
    rate = rate,
    gp_scale = pars$scale,
    gp_shape = pars$shape
  )
}

#' Return levels of an equal-weight mixture of distributions
#'
#' The marginal runoff distribution implied by a distributional learning model
#' is the average of its predictive distributions over the observed peaks (the
#' law of total probability). Its survival function is therefore the average
#' of the component survival functions, which is evaluated here on a grid and
#' inverted directly, without building a mixture object.
#'
#' @param dsts List of distributions (null distributions are ignored).
#' @param return_periods Return periods in years.
#' @param rate Peaks per year.
#' @param n_grid Grid size for the inversion.
#' @returns Numeric vector of return levels (`NA` where `rate * T <= 1`).
#' @export
mixture_return_levels <- function(dsts, return_periods, rate, n_grid = 2000L) {
  keep <- vapply(dsts, function(d) !identical(distionary::pretty_name(d), "Null"), logical(1))
  dsts <- dsts[keep]
  p_exceed <- 1 / (rate * return_periods)
  surv_mix <- function(z) {
    Reduce(`+`, lapply(dsts, function(d) as.numeric(distionary::eval_survival(d, at = z)))) / length(dsts)
  }
  hi <- 1
  target <- min(p_exceed[p_exceed < 1])
  iter <- 0L
  while (surv_mix(hi) > target && iter < 40L) {
    hi <- hi * 2
    iter <- iter + 1L
  }
  zs <- seq(0, hi, length.out = n_grid)
  s <- surv_mix(zs)
  keep <- c(TRUE, diff(s) < 0)
  out <- rep(NA_real_, length(return_periods))
  ok <- p_exceed < 1
  out[ok] <- stats::approx(rev(s[keep]), rev(zs[keep]), xout = p_exceed[ok], rule = 2, ties = "ordered")$y
  out
}

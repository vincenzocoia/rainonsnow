# Copula transport of a cell's marginal, as an alternative to the mixture.
#
# WHY THERE IS AN ALTERNATIVE AT ALL
#
# The cell marginal is currently the equal-weight mixture of that cell's
# peak-hour predictive distributions. Under a shared shape every component has
# the same tail index, so the mixture has that index too -- the CONDITIONAL one.
# The marginal's own index is heavier by the factor CEVI, and the difference
# comes from covariate values beyond the sample, which no amount of learning
# reaches. Simulation puts the shortfall at 20-40% of the 10,000-event return
# level at dependence worth modelling, and it does not shrink with more data
# (see the README section on copula transport).
#
# The copula supplies the missing step. Each peak hour transports its own
# conditional tail back into an estimate of the WHOLE marginal, and the hours
# are then combined as repeated estimates of one quantity rather than mixed.
#
# THE ONE APPROXIMATION, STATED PLAINLY
#
# The transport needs a scalar rank u = F_X(x). The cell has several predictors,
# so there is no such thing directly. What it does have is the forest's own
# one-dimensional summary of each hour: the predictive quantile the tail is
# grafted at, `graft_of`, which is high exactly when the forest thinks that hour
# is primed for a large response. This uses its rank as u.
#
# That is conditioning on the predictive location rather than on the predictor
# vector, and the two agree only if the predictive distributions form a location
# family in it. They do not exactly. The consequence is that the fitted copula
# absorbs whatever the summary loses, so the transported marginal should be read
# alongside the mixture rather than as a drop-in replacement -- which is why it
# is off by default and written to its own file.
#
# WHAT TO CHECK BEFORE BELIEVING IT
#
# Every hour estimates the same marginal, so the spread across hours is a free
# consistency check on the whole construction. `dl_transport_spread()` reports
# it. A cell whose hours disagree by an order of magnitude is telling you the
# conditional model and the copula do not fit together there.

#' Copula parameter from ranks
#'
#' Fits the copula parameter by inverting Kendall's tau, which uses all the
#' cell's hours rather than only its tail and needs no optimiser.
#'
#' @param a,b Numeric vectors of equal length.
#' @param family `"gaussian"`, `"clayton"` or `"survival_clayton"`.
#' @returns A single number: the correlation for `"gaussian"`, theta for the
#'   Clayton families.
#' @examples
#' dl_copula_par(1:20, c(1:15, 20:16), "gaussian")
#' @export
dl_copula_par <- function(a, b, family = "survival_clayton") {
  ok <- is.finite(a) & is.finite(b)
  if (sum(ok) < 5L) {
    return(NA_real_)
  }
  tau <- stats::cor(a[ok], b[ok], method = "kendall")
  if (!is.finite(tau)) {
    return(NA_real_)
  }
  if (identical(family, "gaussian")) {
    return(sin(pi / 2 * tau))
  }
  # Clayton and its survival rotation both have tau = theta / (theta + 2), and
  # both need theta > 0, so a non-positive tau has no Clayton in it.
  if (tau <= 0) {
    return(NA_real_)
  }
  2 * tau / (1 - tau)
}

#' Transport a cell's peak-hour tails into an estimate of its marginal
#'
#' @param df A cell's encoded predictions: one row per peak hour with
#'   `graft_of`, `graft_tail_prob`, `gp_scale` and `gp_shape`, as
#'   [dl_encode_peak_hour_distributions()] writes them.
#' @param observed The response actually observed at each peak hour. Needed to
#'   fit the copula; supply `par` instead to skip it.
#' @param family Copula family; see [cop_h()]. The default is the survival
#'   Clayton, which has upper tail dependence -- the relevant kind when the
#'   predictors and the response are extreme together.
#' @param par Copula parameter. Fitted from `observed` when `NULL`. Unused when
#'   `family` is a function supplying the transform directly.
#' @param location Scalar predictive summary per hour, whose rank becomes `u`.
#'   Defaults to `graft_of`; see the file header for what that assumes.
#' @returns A [marginal_ensemble()], or `NULL` when the cell has too little to
#'   work with.
#' @examples
#' \dontrun{
#' me <- dl_cell_transport_ensemble(cell_rows, observed = cell_rows$runoff_hourly)
#' marginal_ensemble_quantile(me, 1e-4)
#' }
#' @export
dl_cell_transport_ensemble <- function(
  df,
  observed = NULL,
  family = "survival_clayton",
  par = NULL,
  location = NULL
) {
  needed <- c("graft_of", "gp_scale", "gp_shape")
  missing <- setdiff(needed, names(df))
  if (length(missing)) {
    stop(
      "Encoded predictions are missing: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  if (is.null(location)) {
    location <- df$graft_of
  }
  tail_prob <- if ("graft_tail_prob" %in% names(df)) {
    df$graft_tail_prob
  } else {
    rep(0.5, nrow(df))
  }

  # A function standing in for the family supplies the transform itself, so
  # there is no parameter to fit. That is how the simulation lab feeds in a
  # copula that belongs to no family, and how the tests check the identity.
  if (!is.function(family)) {
    if (is.null(par)) {
      if (is.null(observed)) {
        stop(
          "Supply either `observed` (to fit the copula) or `par`.",
          call. = FALSE
        )
      }
      if (length(observed) != nrow(df)) {
        stop("`observed` needs one value per row of `df`.", call. = FALSE)
      }
      par <- dl_copula_par(location, observed, family)
    }
    if (!is.finite(par)) {
      return(NULL)
    }
  }

  keep <- is.finite(location) &
    is.finite(df$graft_of) &
    is.finite(df$gp_scale) &
    df$gp_scale > 0 &
    is.finite(df$gp_shape) &
    is.finite(tail_prob) &
    tail_prob > 0
  if (sum(keep) < 5L) {
    return(NULL)
  }

  # Ranks on (0, 1), open at both ends: the copula's h-function is undefined at
  # 0 and 1, and the largest predictive location should not be treated as the
  # largest one possible.
  u <- rank(location[keep]) / (sum(keep) + 1L)

  out <- tryCatch(
    marginal_ensemble(
      u = u,
      threshold = df$graft_of[keep],
      tail_prob = tail_prob[keep],
      scale = df$gp_scale[keep],
      shape = df$gp_shape[keep],
      family = family,
      par = if (is.function(family)) 0.5 else par
    ),
    error = function(e) NULL
  )
  out
}

#' How far the peak hours disagree about the same marginal
#'
#' Every hour transports to an estimate of the whole marginal, so they should
#' agree. This returns the ratio of the 90th to the 10th percentile of those
#' estimates at one exceedance probability: 1 is perfect agreement, and a large
#' value says the conditional model and the copula are not mutually consistent
#' in that cell.
#'
#' @param me A [marginal_ensemble()].
#' @param p Exceedance probability at which to compare.
#' @returns A single number, or `NA` when too few hours give a usable answer.
#' @examples
#' me <- marginal_ensemble(c(0.3, 0.5, 0.7), c(9, 10, 11), rep(0.3, 3),
#'                         c(2, 2.5, 3), 0.1, par = 0.7)
#' dl_transport_spread(me, 1e-4)
#' @export
dl_transport_spread <- function(me, p = 1e-4) {
  y <- marginal_ensemble_quantile(me, p, "median")
  if (!is.finite(y)) {
    return(NA_real_)
  }
  # Compare the hours on the probability they assign to one level, which needs
  # no per-hour inversion.
  s <- as.vector(marginal_ensemble_paths(me, y))
  s <- s[is.finite(s) & s > 0]
  if (length(s) < 5L) {
    return(NA_real_)
  }
  q <- stats::quantile(s, c(0.1, 0.9), names = FALSE)
  q[2] / q[1]
}

#' Return levels from a transported marginal
#'
#' @param me A [marginal_ensemble()].
#' @param return_periods Numeric vector of return periods, in years.
#' @param events_per_year Peak-over-threshold events per year for the cell.
#' @param combine Combination rule; see [marginal_ensemble_survival()]. The
#'   default pointwise median is the one that does not let a single heavy hour
#'   set the answer.
#' @returns A data frame with `return_period` and `return_level`, `NA` at return
#'   periods below what the transported tail covers.
#' @examples
#' me <- marginal_ensemble(c(0.3, 0.5, 0.7), c(9, 10, 11), rep(0.3, 3),
#'                         c(2, 2.5, 3), 0.1, par = 0.7)
#' dl_transport_return_levels(me, c(10, 100), events_per_year = 3)
#' @export
dl_transport_return_levels <- function(
  me,
  return_periods,
  events_per_year,
  combine = c("median", "mean", "weighted_mean")
) {
  combine <- match.arg(combine)
  p <- 1 / (return_periods * events_per_year)
  data.frame(
    return_period = return_periods,
    return_level = marginal_ensemble_quantile(me, p, combine)
  )
}

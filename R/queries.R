#' Grid of predictor values for querying a fitted analysis
#'
#' Each predictor runs from zero to `mult` times its largest value at the
#' cell's peaks.
#'
#' @param data Training data for one cell.
#' @param predictors Predictor names.
#' @param n Grid points per predictor.
#' @param mult Multiplier on the observed maximum.
#' @returns A tibble with one column per predictor.
#' @export
query_grid <- function(data, predictors, n, mult = 1.1) {
  axes <- lapply(predictors, function(p) {
    seq(0, mult * max(data[[p]], na.rm = TRUE), length.out = n)
  })
  names(axes) <- predictors
  tibble::as_tibble(expand.grid(axes, KEEP.OUT.ATTRS = FALSE))
}

#' Event probability across the predictor grid
#'
#' For each grid point `x` and each event level `z_T`, computes
#' `P(runoff > z_T | x)` from the (tail-treated) predictive distribution.
#' This is the building block for all queries.
#'
#' @param model A fitted `"dstlrn"` object for one cell.
#' @param grid Output of [query_grid()].
#' @param levels Tibble with `return_period` and `return_level`.
#' @param tail The `tail` block of the analysis configuration.
#' @returns The grid crossed with `return_period`, plus `return_level` and
#'   `p_event`.
#' @export
query_event_probability <- function(model, grid, levels, tail) {
  dsts <- dl_apply_tail(predict(model, newdata = grid), tail)
  z <- levels$return_level
  surv <- vapply(dsts, function(d) {
    if (identical(distionary::pretty_name(d), "Null")) {
      return(rep(NA_real_, length(z)))
    }
    as.numeric(distionary::eval_survival(d, at = z))
  }, numeric(length(z)))
  surv <- matrix(surv, nrow = length(z))
  out <- grid[rep(seq_len(nrow(grid)), each = length(z)), ]
  out$return_period <- rep(levels$return_period, times = nrow(grid))
  out$return_level <- rep(z, times = nrow(grid))
  out$p_event <- as.numeric(surv)
  out
}

#' Trigger thresholds: how much of the target predictor sets off the event
#'
#' For each event level and each combination of the other predictors, the
#' smallest value of `target` at which `P(runoff > z_T | x)` reaches `prob`,
#' linearly interpolated between grid points. `NA` means the probability is
#' not reached anywhere on the grid (i.e. not within the observed range).
#'
#' @param event_tbl Output of [query_event_probability()].
#' @param target The predictor being solved for (e.g. rainfall).
#' @param probs Probability levels.
#' @returns A tibble with the conditioning predictors, `return_period`,
#'   `return_level`, `prob` and `trigger`.
#' @export
query_trigger <- function(event_tbl, target, probs) {
  others <- setdiff(
    names(event_tbl),
    c(target, "return_period", "return_level", "p_event")
  )
  event_tbl |>
    dplyr::group_by(dplyr::across(dplyr::all_of(c(others, "return_period", "return_level")))) |>
    dplyr::group_modify(function(d, key) {
      d <- d[order(d[[target]]), ]
      tibble::tibble(
        prob = probs,
        trigger = vapply(probs, function(p) {
          first_crossing(d[[target]], d$p_event, p)
        }, numeric(1))
      )
    }) |>
    dplyr::ungroup()
}

# Smallest x at which y first reaches level, interpolating linearly between
# the last point below and the first point at or above.
first_crossing <- function(x, y, level) {
  y[is.na(y)] <- -Inf
  i <- which(y >= level)[1]
  if (is.na(i)) return(NA_real_)
  if (i == 1) return(x[1])
  x0 <- x[i - 1]; x1 <- x[i]; y0 <- y[i - 1]; y1 <- y[i]
  if (!is.finite(y0)) return(x1)
  x0 + (level - y0) / (y1 - y0) * (x1 - x0)
}

#' Likeliest drivers of a T-year event
#'
#' The density of the predictors given that the event occurred,
#' `f(x | runoff > z_T)`, proportional to `P(runoff > z_T | x) f(x)`, is
#' normalised numerically over the grid. With two predictors it also returns
#' the density of the target conditional on the other predictor as well,
#' `f(target | runoff > z_T, other)`, normalised over the target for each
#' value of the other predictor.
#'
#' @param event_tbl Output of [query_event_probability()].
#' @param drivers A `"drivers"` object from [fit_drivers()].
#' @param target Target predictor.
#' @returns `event_tbl` with `f_prior` (drivers density), `f_given_event`, and
#'   for two predictors `f_prior_cond` and `f_given_event_cond`.
#' @export
query_likeliest <- function(event_tbl, drivers, target) {
  preds <- drivers$predictors
  grid <- dplyr::distinct(event_tbl[preds])
  grid$f_prior <- eval_drivers_density(drivers, grid)
  if (length(preds) == 2) {
    other <- setdiff(preds, target)
    grid$f_prior_cond <- eval_drivers_density(drivers, grid, given = other)
  }
  cell_area <- prod(vapply(preds, function(p) {
    v <- sort(unique(grid[[p]]))
    if (length(v) > 1) v[2] - v[1] else 1
  }, numeric(1)))
  target_step <- diff(sort(unique(grid[[target]])))[1]
  out <- dplyr::left_join(event_tbl, grid, by = preds) |>
    dplyr::group_by(.data$return_period) |>
    dplyr::mutate(
      f_given_event = normalise_density(.data$p_event * .data$f_prior, cell_area)
    )
  if (length(preds) == 2) {
    out <- out |>
      dplyr::group_by(.data$return_period, .data[[other]]) |>
      dplyr::mutate(
        f_given_event_cond = normalise_density(
          .data$p_event * .data$f_prior_cond, target_step
        )
      )
  }
  dplyr::ungroup(out)
}

normalise_density <- function(v, step) {
  v[!is.finite(v)] <- 0
  total <- sum(v) * step
  if (total <= 0) return(rep(NA_real_, length(v)))
  v / total
}

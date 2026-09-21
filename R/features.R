#' Predictor registry
#'
#' Every predictor an analysis can use is defined here, as a function of one
#' cell's hourly table (sorted by `date`) that returns a numeric vector with one
#' value per hour. Defining predictors on the hourly series, rather than on the
#' peaks, is what allows ones that need history: accumulations, rolling sums,
#' or a snowpack water balance.
#'
#' To add a predictor, add an entry here and name it under `predictors:` in an
#' `analysis.yaml`. Units follow `scripts/2-tablify_spatial_eo.r` (mm/h for
#' rates; ERA5-Land depths are in metres).
#'
#' @returns A named list of functions.
#' @examples
#' names(feature_registry())
#' @export
feature_registry <- function() {
  list(
    rainfall_hourly = function(h) h$rainfall_hourly,
    snowmelt_hourly = function(h) h$snowmelt_hourly,
    # Snow water equivalent in mm: a first proxy for water stored in the pack.
    swe_mm = function(h) h$snow_depth_water_equivalent * 1000,
    # Rain accumulated over the preceding day, including the current hour.
    rainfall_24h = function(h) rolling_sum(h$rainfall_hourly, 24L),
    snowmelt_24h = function(h) rolling_sum(h$snowmelt_hourly, 24L)
  )
}

#' Trailing rolling sum
#'
#' @param x Numeric vector in time order.
#' @param k Window length (number of steps, including the current one).
#' @returns Numeric vector the same length as `x`; the first `k - 1` values sum
#'   over the shorter available window.
#' @examples
#' rolling_sum(1:5, 2)
#' @export
rolling_sum <- function(x, k) {
  x[is.na(x)] <- 0
  cs <- cumsum(x)
  lagged <- c(rep(0, k), cs)[seq_along(cs)]
  cs - lagged
}

#' Build an analysis training table
#'
#' Computes the requested predictors on each cell's hourly series, then keeps
#' the POT peak hours.
#'
#' @param hourly Hourly table from `scripts/2-tablify_spatial_eo.r`.
#' @param peaks POT peaks from `scripts/3-pot_spatial_eo.r`.
#' @param predictors Names from [feature_registry()].
#' @param response Response column.
#' @returns A tibble with `cell_id`, `x`, `y`, `date`, the response and the
#'   predictors, one row per peak.
#' @export
build_training <- function(hourly, peaks, predictors, response = "runoff_hourly") {
  reg <- feature_registry()[predictors]
  keys <- c("cell_id", "x", "y", "date")
  feats <- hourly |>
    dplyr::arrange(.data$cell_id, .data$date) |>
    dplyr::group_by(.data$cell_id) |>
    dplyr::group_modify(function(h, key) {
      vals <- lapply(reg, function(f) f(h))
      tibble::as_tibble(c(list(date = h$date), vals))
    }) |>
    dplyr::ungroup()
  out <- peaks[c(keys, response)]
  out <- dplyr::left_join(out, feats, by = c("cell_id", "date"))
  out[c(keys, response, predictors)]
}

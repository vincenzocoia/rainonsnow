#' Tabulate marginal return levels on a POT event-frequency axis
#'
#' Calls `distionary::enframe_return()` at abscissae `return_periods *
#' num_events_per_year`, i.e. return periods expressed in **events** (per the
#' PoT rate), not calendar years. Adds columns `return_period_years` (the same
#' grid as `return_periods`) and `num_events_per_year` so downstream code can
#' filter or plot using reporting years while `return_period` stays event-scaled.
#'
#' @param marginal A distribution object passed to `distionary::enframe_return()`
#'   as its first argument.
#' @param num_events_per_year Estimated PoT rate for the cell (scalar):
#'   expected peaks per calendar year.
#' @param return_periods Reporting return periods in **years**, in order; defaults
#'   to [rp_reporting()]. Values are multiplied by `num_events_per_year` before
#'   passing to `enframe_return()`, which names that column `return_period`
#'   (event units).
#'
#' @returns A tibble: columns from `enframe_return()` (including `return_period`
#'   in event units and `return`), plus `return_period_years` and
#'   `num_events_per_year`.
#'
#' @seealso [rp_reporting()]
#' @export
#' @examples
#' \dontrun{
#' enframe_at_events(my_marginal, num_events_per_year = 4.1)
#' enframe_at_events(my_marginal, 4.1, return_periods = c(2, 10, 100))
#' }
enframe_at_events <- function(
    marginal,
    num_events_per_year,
    return_periods = rp_reporting()
) {
  nep <- num_events_per_year
  distionary::enframe_return(
    marginal,
    at = return_periods * nep,
    arg_name = "return_period"
  ) |>
    dplyr::mutate(
      return_period_years = return_periods,
      num_events_per_year = nep
    )
}

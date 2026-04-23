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

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
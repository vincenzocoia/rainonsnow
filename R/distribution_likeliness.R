#' Likeliness of a distribution to produce a specified value
#'
#' For a set of distributions, calculate the likeliness of each
#' distribution to produce a specified value.
#' @param distributions A list of distributions.
#' @param value The value to calculate the likeliness for.
#' @returns A vector of likeliness values, one for each distribution.
#' @details
#' The likeliness is calculated as the density of the distribution at the
#' specified value divided by the sum of the densities of all distributions.
#'
#' If a density evaluates to `NA`, it is ignored in the sum.
#' @examples
#' distributions <- list(
#'   distionary::dst_normal(mean = 0, sd = 1),
#'   distionary::dst_normal(mean = 1, sd = 1),
#'   distionary::dst_null()
#' )
#' distribution_likeliness(distributions, 0)
#' @export
distribution_likeliness <- function(distributions, value) {
  densities <- purrr::map_dbl(
    distributions,
    distionary::eval_density,
    at = value
  )
  densities / sum(densities, na.rm = TRUE)
}

#' Regular grid from a scatter
#'
#' For a scattering of x-y points, this function makes a regular x-y grid
#' based on the range of each variable.
#' @param x,y Names of the x and y variables (unquoted); they will be evaluated
#'   in the `data` data mask.
#' @param data A possible data frame to evaluate the x and y variables in; NULL
#'   if not needed (i.e., grab vectors from the enclosing environment).
#' @param size The number of grid points in each dimension; vector of length 1
#'   for the same size in both dimensions, or vector of length 2 corresponding
#'   to (x, y).
#' @param mult A multiplier for the range of each variable; vector of length 1
#'   for the same size in both dimensions, or vector of length 2 corresponding
#'   to (x, y). NOTE: expansion only happens in the positive direction (bottom
#'   left corner stays the same) because the focus is on extremes.
#' @returns A data frame with columns `x` and `y`.
#' @examples
#' set.seed(123)
#' x <- rnorm(100)
#' y <- rnorm(100)
#' grid <- grid_from_scatter(x, y, size = 10, mult = 1.1)
#' plot(grid$x, grid$y)
#' points(x, y, col = "red", pch = 19)
#' legend(
#'   "bottomright",
#'   legend = c("Grid", "Scatter"),
#'   col = c("black", "red"),
#'   pch = c(1, 19)
#' )
#' @export
grid_from_scatter <- function(x, y, data = NULL, size = 30, mult = 1) {
  checkmate::assert_numeric(size, any.missing = FALSE)
  checkmate::assert_numeric(mult, any.missing = FALSE)
  size <- vctrs::vec_recycle(size, size = 2)
  mult <- vctrs::vec_recycle(mult, size = 2)
  xquo <- rlang::enquo(x)
  yquo <- rlang::enquo(y)
  x <- rlang::eval_tidy(xquo, data = data)
  y <- rlang::eval_tidy(yquo, data = data)
  x_rng <- range(x, na.rm = TRUE)
  y_rng <- range(y, na.rm = TRUE)
  x_rng <- x_rng[1] + c(0, 1) * diff(x_rng) * mult[1]
  y_rng <- y_rng[1] + c(0, 1) * diff(y_rng) * mult[2]
  tidyr::expand_grid(
    x = seq(x_rng[1], x_rng[2], length.out = size[1]),
    y = seq(y_rng[1], y_rng[2], length.out = size[2])
  )
}
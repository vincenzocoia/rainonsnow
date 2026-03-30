#' Predict from a distributional learning model
#'
#' `predict()` is an S3 generic from the stats package. Methods are provided
#' for distributional learning objects inheriting from class `"dstlrn"`.
#'
#' @param object A distributional learning object.
#' @param newdata A data frame of new data to predict. If `NULL`, the training
#'   data is used.
#' @param ... Additional arguments passed on to specific methods.
#' @returns A list of probability distributions, one for each row of `newdata`.
#' @rdname predict.dstlrn
#' @importFrom stats predict
#' @export
predict.dstlrn <- function(object, newdata = NULL, ...) {
  rlang::abort(
    "Don't know how to make predictions from this distributional learning object."
  )
}

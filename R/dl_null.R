#' Make a null distributional learning object
#' 
#' Useful for when a distributional learning method fails
#' 
#' @param data A data frame of training data. Useful if you want
#'   to run `predict()` without specifying `newdata`,
#'   in which case the number of rows of the data frame will be used.
#' @returns A distributional learning object with subclass
#'   "dl_null".
#' @examples
#' dl_null()
#' @export
dl_null <- function(data = NULL) {
  new_dstlrn(list(training = data), subclass = "dl_null")
}

#' @describeIn predict.dstlrn Predict from a null distributional learning model.
#' @export
predict.dl_null <- function(object, newdata, ...) {
  n <- nrow(newdata)
  rep(distionary::dst_null(), n)
}
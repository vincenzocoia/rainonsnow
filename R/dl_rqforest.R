#' Make a quantile regression forest distributional learning model
#' 
#' @param data A data frame.
#' @param yname The name of the response variable.
#' @param xnames The names of the predictor variables.
#' @param ... Other arguments passed down to `quantregForest::quantregForest()`.
#' @param na_action The action to take if there are missing values in the data;
#'   either "drop" to drop rows with missing values, or "null" to return a
#'   "dl_null" object.
#' @param min_obs The minimum number of observations required to fit the model.
#'   If the number of observations is less than `min_obs`, a "dl_null" object
#'   is returned.
#' @returns A distributional learning model with subclass "dl_rqforest".
#' @examples
#' set.seed(123)
#' df <- tibble::tibble(x = rnorm(100), y = x + rnorm(100))
#' dl_rqforest(data = df, yname = "y", xnames = "x")
#' @export
dl_rqforest <- function(data, yname, xnames, ..., na_action = c("drop", "null"), min_obs = 5) {
  na_action <- rlang::arg_match(na_action)
  data <- data[append(xnames, yname)]
  data2 <- tidyr::drop_na(data)
  if (na_action == "null" && nrow(data2) < nrow(data)) {
    return(dl_null())
  }
  data <- data2
  if (nrow(data) < min_obs) {
    return(dl_null(data))
  }
  y <- data[[yname]]
  x <- data[xnames]
  base_model <- try(
    quantregForest::quantregForest(x = x, y = y, ...),
    silent = TRUE
  )
  if (inherits(base_model, "try-error")) {
    return(dl_null(data))
  }
  res <- list(
    base = base_model,
    yname = yname,
    xnames = xnames,
    training = data
  )
  new_dstlrn(res, subclass = "dl_rqforest")
}

#' @describeIn predict.dstlrn Predict from a quantile regression forest
#'   distributional learning model.
#' @param ... Other arguments passed down to `predict.quantregForest()`, except
#'   for the `what` argument which is used internally to create probability
#'   distributions.
#' @export
predict.dl_rqforest <- function(object, newdata = NULL, ...) {
  base_model <- object[["base"]]
  if (is.null(newdata)) {
    newdata <- object[["training"]]
  }
  newdata <- newdata[object[["xnames"]]]
  n <- nrow(newdata)
  lgl_na <- df_rows_have_missing(newdata)
  res <- rep(list(distionary::dst_null()), n)
  if (all(lgl_na)) {
    return(res)
  }
  non_nulls <- predict(
    base_model,
    newdata = newdata[!lgl_na, ],
    what = stats::ecdf,
    ...
  )
  res[!lgl_na] <- lapply(non_nulls, ecdf_to_dst)
  res
}
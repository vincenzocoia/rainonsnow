#' Print a distributional learning object
#'
#' @param x A distributional learning object.
#' @param ... Unused.
#' @returns `x`, invisibly.
#' @export
print.dstlrn <- function(x, ...) {
  subclass <- setdiff(class(x), "dstlrn")[1]
  training <- x[["training"]]
  if (is.na(subclass)) {
    subclass <- "dstlrn"
  }

  cat("<", subclass, ">\n", sep = "")

  if (!is.null(x[["yname"]])) {
    cat("response: ", x[["yname"]], "\n", sep = "")
  }

  if (!is.null(x[["xnames"]])) {
    cat("predictors: ", length(x[["xnames"]]), "\n", sep = "")
  }

  if (is.data.frame(training)) {
    cat("training rows: ", nrow(training), "\n", sep = "")
  }

  if (identical(subclass, "dl_null")) {
    cat("status: null model\n", sep = "")
  }

  invisible(x)
}

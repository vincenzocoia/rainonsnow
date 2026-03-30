#' Make a distributional learning object
#' 
#' @param l A list to use as the basis for the object.
#' @param ... Other arguments to add as attributes to the object.
#' @param subclass The subclass of the distributional learning object.
#' @returns A distributional learning object with subclass `subclass`.
#' @examples
#' new_dstlrn(list())
#' new_dstlrn(list(a = 1, b = 2), subclass = "my_subclass")
#' @export
new_dstlrn <- function(x, ..., subclass = character()) {
  structure(x, ..., class = c(subclass, "dstlrn"))
}
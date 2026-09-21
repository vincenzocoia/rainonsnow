#' Distributional learning model registry
#'
#' Maps the `model: type:` field of an `analysis.yaml` to a fitting function
#' with signature `function(data, yname, xnames, ...)` that returns a
#' `"dstlrn"` object with a `predict()` method. To add a model (e.g. a wrapper
#' that transports another model's predictions through a copula), write the
#' constructor and `predict()` method, then register it here.
#'
#' Each entry also records the number of predictors the model supports.
#'
#' @returns A named list; each element has `fit` (function) and `max_p`
#'   (integer).
#' @export
dl_model_registry <- function() {
  list(
    llqr = list(fit = dl_llqr, max_p = 1L),
    rqforest = list(fit = dl_rqforest, max_p = Inf)
  )
}

#' Fit the configured distributional learning model to one cell
#'
#' @param data Training data for one cell.
#' @param cfg Analysis configuration from [read_analysis()].
#' @returns A `"dstlrn"` object.
#' @export
dl_fit <- function(data, cfg) {
  entry <- dl_model_registry()[[cfg$model$type]]
  if (length(cfg$predictors) > entry$max_p) {
    rlang::abort(paste0(
      "Model '", cfg$model$type, "' supports at most ", entry$max_p,
      " predictor(s); analysis '", cfg$name, "' has ",
      length(cfg$predictors), "."
    ))
  }
  args <- c(
    list(data = data, yname = cfg$response, xnames = cfg$predictors),
    cfg$model$args
  )
  do.call(entry$fit, args)
}

#' Apply the configured tail treatment to predictive distributions
#'
#' With `tail: type: gp`, a generalized Pareto tail is grafted onto each step
#' distribution (see [fit_and_graft_gp()]), so that exceedance probabilities
#' beyond the largest training value are not forced to zero. Null or
#' non-finite distributions, and grafts that fail, are passed through as null
#' distributions.
#'
#' @param dsts List of predictive distributions.
#' @param tail The `tail` block of the configuration.
#' @returns A list of distributions the same length as `dsts`.
#' @export
dl_apply_tail <- function(dsts, tail) {
  if (identical(tail$type, "none")) {
    return(dsts)
  }
  lapply(dsts, function(d) {
    if (!identical(distionary::pretty_name(d), "Finite")) {
      return(distionary::dst_null())
    }
    tryCatch(
      suppressWarnings(fit_and_graft_gp(d, adaptive_threshold = tail$adaptive_threshold)),
      error = function(e) distionary::dst_null()
    )
  })
}

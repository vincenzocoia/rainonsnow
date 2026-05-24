#' Compact storage for grafted GP-tail predictions
#'
#' @description
#' Grafted distributions from [fit_and_graft_gp()] are large when stored as a
#' list column. This module stores only the GPD tail (`graft_of`, `gp_scale`,
#' `gp_shape`) and rebuilds `distribution_gp` from `distribution_forest` when
#' needed.
#'
#' @section Which function to use:
#'
#' * **Save or load the predictions RDS** (scripts, apps): [dl_write_peak_hour_predictions()]
#'   and [dl_read_peak_hour_predictions()]. The file on disk is always the
#'   compact form (no `distribution_gp` column).
#' * **Already in memory, about to write**: [dl_encode_peak_hour_distributions()]
#'   after fitting; diagnostics can run on the full tibble *before* encoding.
#' * **Already loaded a compact tibble**: [dl_decode_peak_hour_distributions()]
#'   before anything that needs `distribution_gp` (e.g. [mix2()], [dl_pp_long()]).
#'
#' @name dl_graft_codec
NULL

graft_gp_parameters <- function(dst) {
  checkmate::assert_class(dst, "dst")
  pars <- distionary::parameters(dst)
  dists <- pars$distributions
  checkmate::assert_list(dists, len = 2L, types = "dst")

  tail_pars <- distionary::parameters(dists[[2L]])
  gp_pars <- distionary::parameters(tail_pars$distribution)
  ss <- unlist(gp_pars[c("scale", "shape")])

  list(
    graft_of = tail_pars$shift,
    gp_scale = unname(ss[["scale"]]),
    gp_shape = unname(ss[["shape"]])
  )
}

reconstruct_graft_gp <- function(
  distribution_forest,
  graft_of,
  gp_scale,
  gp_shape,
  include = FALSE
) {
  checkmate::assert_class(distribution_forest, "dst")
  checkmate::assert_number(graft_of)
  checkmate::assert_number(gp_scale, lower = 0)
  checkmate::assert_number(gp_shape)

  tail <- distionary::dst_gp(gp_scale, gp_shape) + graft_of
  distplyr::graft_right(
    distribution_forest,
    tail,
    of = graft_of,
    include = include
  )
}

#' Encode peak-hour distributional-learning predictions (in memory)
#'
#' Drops `distribution_gp` and adds `graft_of`, `gp_scale`, and `gp_shape`. Prefer
#' [dl_write_peak_hour_predictions()] when writing to disk.
#'
#' @param peak_hour_distributions Tibble with `distribution_gp` (and
#'   `distribution_forest`).
#' @return Compact tibble without `distribution_gp`.
#' @seealso [dl_decode_peak_hour_distributions()], [dl_write_peak_hour_predictions()]
#' @export
dl_encode_peak_hour_distributions <- function(peak_hour_distributions) {
  checkmate::assert_data_frame(peak_hour_distributions)
  checkmate::assert_names(
    names(peak_hour_distributions),
    must.include = "distribution_gp"
  )
  gp_enc <- purrr::map_dfr(
    peak_hour_distributions$distribution_gp,
    graft_gp_parameters
  )
  peak_hour_distributions |>
    dplyr::select(-"distribution_gp") |>
    dplyr::bind_cols(gp_enc)
}

#' Decode compact peak-hour predictions (in memory)
#'
#' Rebuilds `distribution_gp` for downstream use ([mix2()], [dl_pp_long()], etc.).
#' Prefer [dl_read_peak_hour_predictions()] when reading from disk.
#'
#' @param encoded Tibble from [dl_encode_peak_hour_distributions()] or
#'   [dl_read_peak_hour_predictions()] (before decoding).
#' @param include Passed to `graft_right()`; default `FALSE` (same as
#'   [fit_and_graft_gp()]).
#' @return Tibble with `distribution_gp`.
#' @seealso [dl_encode_peak_hour_distributions()], [dl_read_peak_hour_predictions()]
#' @export
dl_decode_peak_hour_distributions <- function(encoded, include = FALSE) {
  checkmate::assert_data_frame(encoded)
  checkmate::assert_names(
    names(encoded),
    must.include = c(
      "distribution_forest",
      "graft_of",
      "gp_scale",
      "gp_shape"
    )
  )
  encoded |>
    dplyr::mutate(
      distribution_gp = purrr::pmap(
        list(
          .data$distribution_forest,
          .data$graft_of,
          .data$gp_scale,
          .data$gp_shape
        ),
        reconstruct_graft_gp,
        include = include
      )
    )
}

#' @rdname dl_encode_peak_hour_distributions
#' @export
dl_predictions_is_encoded <- function(x) {
  is.data.frame(x) &&
    "graft_of" %in% names(x) &&
    !"distribution_gp" %in% names(x)
}

#' Write compact peak-hour predictions to disk
#'
#' Encodes `distribution_gp` to numeric columns if present, then writes RDS.
#'
#' @param peak_hour_distributions Full or already-encoded predictions tibble.
#' @param path Path for `saveRDS()`.
#' @param ... Passed to [saveRDS()].
#' @return Invisibly, the encoded tibble that was written.
#' @seealso [dl_read_peak_hour_predictions()]
#' @export
dl_write_peak_hour_predictions <- function(
  peak_hour_distributions,
  path,
  ...
) {
  encoded <- if (dl_predictions_is_encoded(peak_hour_distributions)) {
    peak_hour_distributions
  } else {
    dl_encode_peak_hour_distributions(peak_hour_distributions)
  }
  saveRDS(encoded, file = path, ...)
  invisible(encoded)
}

#' Read peak-hour predictions and restore `distribution_gp`
#'
#' Reads the compact RDS written by [dl_write_peak_hour_predictions()] and
#' rebuilds `distribution_gp`. If the file already contains `distribution_gp`
#' (legacy format), it is returned unchanged.
#'
#' @param path Path passed to [readRDS()].
#' @param decode If `TRUE` (default), rebuild `distribution_gp` when the file is
#'   encoded.
#' @param ... Passed to [readRDS()].
#' @return Predictions tibble with `distribution_gp` when `decode = TRUE`.
#' @seealso [dl_write_peak_hour_predictions()]
#' @export
dl_read_peak_hour_predictions <- function(path, decode = TRUE, ...) {
  x <- readRDS(file = path, ...)
  if (decode && dl_predictions_is_encoded(x)) {
    x <- dl_decode_peak_hour_distributions(x)
  }
  x
}

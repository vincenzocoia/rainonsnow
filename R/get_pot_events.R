#' Get Peaks Over Threshold Events
#'
#' Identify events where the flow exceeds a given threshold.
#'
#' @param data A data frame containing the flow and date columns.
#' @param threshold The threshold value; either a function that takes a vector
#'   of flow values and returns a single threshold value, or a single numeric
#'   threshold value.
#' @param ... Additional arguments, currently unused.
#' @param flow_col The name of the flow column (character).
#' @param date_col The name of the date column (character).
#' @param min_gap The smallest number of observations in between two consecutive
#'   "runs" (i.e., contiguous threshold exceedances) to be considered a single
#'   event. Non-negative numeric. If 1, no merging is done, because each run
#'   will always have at least one non-exceedance between them.
#'
#' @returns The original data frame (as a tibble), filtered to the highest
#'   flow value for each event. If the max flow is sustained for more than one
#'   row, the first row is returned.
#'
#' The threshold value is also available, stored as an attribute of the returned
#'   data frame, named "threshold".
#' @examples
#' data <- tibble::tibble(
#'   date = as.Date(paste0("2020-01-", 1:30)),
#'   flow = c(
#'     1, 2, 3, 4, 5,
#'     5, 3, 2, 1, 0,
#'     1, 2, 5, 4, 4,
#'     5, 3, 2, 1, 0,
#'     1, 2, 3, 4, 5,
#'     4, 3, 2, 1, 0
#'   )
#' )
#' plot(data$flow, type = "l")
#' get_pot_events(data, threshold = 2)
#'
#' ## merge runs separated by fewer than three observations:
#' x <- tibble::tibble(
#'   date = lubridate::ymd_hms("2020-01-01 00:00:00") +
#'     lubridate::dhours(6) * 0:12,
#'   flow = c(4, 4, 0, 5, 5, 0, 0, 2, 0, 0, 0, 6, 6)
#' )
#' ggplot2::ggplot(x, ggplot2::aes(date, flow)) +
#'   ggplot2::geom_line() +
#'   ggplot2::geom_point()
#' get_pot_events(x, threshold = 1, min_gap = 1)
#' get_pot_events(x, threshold = 1, min_gap = 2)
#' get_pot_events(x, threshold = 1, min_gap = 3)
#' get_pot_events(x, threshold = 1, min_gap = 4)
#' @export
get_pot_events <- function(data,
                           threshold,
                           ...,
                           flow_col = "flow",
                           date_col = "date",
                           min_gap = 0) {
  rlang::check_dots_empty()
  checkmate::assert_data_frame(data)
  checkmate::assert_character(flow_col)
  checkmate::assert_character(date_col)
  checkmate::assert_true(all(c(flow_col, date_col) %in% names(data)))
  date_vec <- data[[date_col]]
  checkmate::assert_true(inherits(date_vec, "Date") || inherits(date_vec, "POSIXt"))
  checkmate::assert_numeric(data[[flow_col]])
  checkmate::assert_number(min_gap, lower = 0, finite = TRUE)

  thr <- if (is.function(threshold)) {
    threshold(data[[flow_col]])
  } else {
    threshold
  }
  checkmate::assert_number(thr, na.ok = TRUE)
  if (is.na(thr)) {
    df <- dplyr::slice_head(tibble::as_tibble(data), n = 0)
    attr(df, "threshold") <- thr
    return(df)
  }

  prep <- tibble::as_tibble(data) |>
    dplyr::arrange(.data[[date_col]]) |>
    dplyr::mutate(
      .idx = dplyr::row_number(),
      .exceed = !is.na(.data[[flow_col]]) & .data[[flow_col]] > thr,
      .pot_event = dplyr::consecutive_id(.exceed)
    ) |>
    dplyr::filter(.exceed)

  drop_prep_meta <- function(x) {
    dplyr::select(x, -dplyr::any_of(c(".pot_event", ".exceed", ".idx")))
  }

  if (nrow(prep) == 0L) {
    df <- drop_prep_meta(prep)
    attr(df, "threshold") <- thr
    return(df)
  }

  runs <- dplyr::summarise(
    prep,
    idx_first = min(.idx),
    idx_last = max(.idx),
    .by = .pot_event
  ) |>
    dplyr::arrange(idx_first)

  merged_key <- if (min_gap <= 0 || nrow(runs) < 2L) {
    runs$.pot_event
  } else {
    .merge_pot_run_keys(runs, min_gap)
  }

  runs <- dplyr::mutate(runs, .merged_event = merged_key)

  df <- dplyr::left_join(
    prep,
    dplyr::select(runs, ".pot_event", ".merged_event"),
    by = ".pot_event"
  ) |>
    dplyr::group_by(.merged_event) |>
    dplyr::arrange(
      dplyr::desc(.data[[flow_col]]),
      .data[[date_col]]
    ) |>
    dplyr::slice_head(n = 1) |>
    dplyr::ungroup() |>
    dplyr::select(-dplyr::any_of(c(".pot_event", ".merged_event", ".exceed", ".idx")))

  attr(df, "threshold") <- thr
  df
}

.merge_pot_run_keys <- function(runs, min_gap) {
  n <- nrow(runs)
  keys <- integer(n)
  keys[[1L]] <- 1L
  cluster <- 1L
  cur_last <- runs$idx_last[[1L]]
  for (i in seq.int(2L, n)) {
    gap_obs <- runs$idx_first[[i]] - cur_last - 1L
    if (gap_obs < min_gap) {
      keys[[i]] <- cluster
      cur_last <- max(cur_last, runs$idx_last[[i]])
    } else {
      cluster <- cluster + 1L
      keys[[i]] <- cluster
      cur_last <- runs$idx_last[[i]]
    }
  }
  keys
}

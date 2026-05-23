#' P–P calibration table for distributional-learning predictions
#'
#' @param peak_hour_distributions Tibble with `cell_id`, `x`, `y`,
#'   `runoff_hourly`, `distribution_forest`, and `distribution_gp`.
#' @return Long tibble with `p_empirical`, `p_model`, and `model` (`forest` / `gp`).
#' @export
dl_pp_long <- function(peak_hour_distributions) {
  peak_hour_distributions |>
    dplyr::group_by(cell_id, x, y) |>
    dplyr::mutate(
      p_model_forest = purrr::map2_dbl(distribution_forest, runoff_hourly, eval_cdf),
      p_empirical_forest = uscore(p_model_forest),
      p_model_gp = purrr::map2_dbl(distribution_gp, runoff_hourly, eval_cdf),
      p_empirical_gp = uscore(p_model_gp)
    ) |>
    dplyr::ungroup() |>
    dplyr::select(cell_id, x, y, dplyr::starts_with("p_")) |>
    tidyr::pivot_longer(
      dplyr::starts_with("p_"),
      names_to = c(".value", "model"),
      names_pattern = "(p_.*)_(.*)"
    ) |>
    dplyr::arrange(cell_id, model, p_empirical)
}

#' Quantile skill scores vs cell-wise empirical marginal
#'
#' @param peak_hour_distributions Same as [dl_pp_long()].
#' @param peaks_dat POT peaks with `cell_id`, `x`, `y`, and `runoff_hourly`
#'   (all hours per cell, for the null marginal).
#' @return Tibble with `cell_id`, `x`, `y`, `tau`, `qscore_null`,
#'   `qscore_model`, and `skill_score`.
#' @export
dl_skill_scores <- function(peak_hour_distributions, peaks_dat) {
  qscores_model <- peak_hour_distributions |>
    dplyr::mutate(
      df = purrr::map(
        distribution_forest,
        enframe_quantile,
        at = 1:99 / 100,
        arg_name = "tau"
      )
    ) |>
    tidyr::unnest(df) |>
    dplyr::group_by(cell_id, x, y, tau) |>
    dplyr::summarise(
      qscore_model = mean(quantile_score(
        runoff_hourly,
        xhat = quantile,
        tau = tau
      )),
      .groups = "drop"
    )

  null_model <- peaks_dat |>
    dplyr::group_by(cell_id, x, y) |>
    dplyr::summarise(runoff_hourly = list(runoff_hourly), .groups = "drop") |>
    dplyr::mutate(marginal = purrr::map(runoff_hourly, dst_empirical))

  null_quantiles <- null_model |>
    dplyr::mutate(
      df = purrr::map(marginal, enframe_quantile, at = 1:99 / 100, arg_name = "tau")
    ) |>
    dplyr::select(!marginal) |>
    tidyr::unnest(df)

  qscores_null <- null_quantiles |>
    dplyr::mutate(
      qscore_null = purrr::pmap_dbl(
        list(runoff_hourly, quantile, tau),
        \(y, q, p) mean(quantile_score(y, xhat = q, tau = p))
      )
    ) |>
    dplyr::select(x, y, tau, qscore_null)

  dplyr::left_join(qscores_null, qscores_model, by = c("x", "y", "tau")) |>
    dplyr::mutate(skill_score = 1 - qscore_model / qscore_null)
}

#' Median quantile skill per cell
#'
#' @param skill_tbl Output of [dl_skill_scores()].
#' @return Tibble with `cell_id`, `x`, `y`, and `skill_median`.
#' @export
dl_map_summary <- function(skill_tbl) {
  skill_tbl |>
    dplyr::group_by(cell_id, x, y) |>
    dplyr::summarise(
      skill_median = median(skill_score, na.rm = TRUE),
      .groups = "drop"
    )
}

#' Build full-grid distributional-learning diagnostics
#'
#' @param peak_hour_distributions Tibble passed to [dl_pp_long()].
#' @param peaks_dat POT peaks passed to [dl_skill_scores()].
#' @return List with `pp`, `skill`, and `map_summary` tibbles.
#' @export
dl_build_diagnostics <- function(peak_hour_distributions, peaks_dat) {
  skill <- dl_skill_scores(peak_hour_distributions, peaks_dat)
  list(
    pp = dl_pp_long(peak_hour_distributions),
    skill = skill,
    map_summary = dl_map_summary(skill)
  )
}

#' Default diagnostics RDS basename (under `derived/`)
#'
#' @export
dl_diagnostics_basename <- function() {
  "era5_land_hourly_alps_dl_diagnostics.rds"
}

#' Write distributional-learning diagnostics to RDS
#'
#' @param diagnostics List from [dl_build_diagnostics()].
#' @param path Output path.
#' @return Invisibly returns `path`.
#' @export
dl_write_diagnostics <- function(diagnostics, path) {
  if (!all(c("pp", "skill", "map_summary") %in% names(diagnostics))) {
    stop("diagnostics must contain pp, skill, and map_summary.", call. = FALSE)
  }
  saveRDS(diagnostics, path)
  invisible(path)
}

#' Read diagnostics RDS written by script 4
#'
#' @param path Path to diagnostics RDS.
#' @return List with `pp`, `skill`, and `map_summary`, or an object of class
#'   `dl_diagnostics_read_error`.
#' @export
dl_read_diagnostics <- function(path) {
  if (!file.exists(path)) {
    return(structure(
      list(message = paste("File not found:", path)),
      class = "dl_diagnostics_read_error"
    ))
  }
  out <- tryCatch(readRDS(path), error = function(e) e)
  if (inherits(out, "error")) {
    return(structure(
      list(message = conditionMessage(out)),
      class = "dl_diagnostics_read_error"
    ))
  }
  if (!is.list(out) || !all(c("pp", "skill", "map_summary") %in% names(out))) {
    return(structure(
      list(message = "Invalid diagnostics RDS (expected pp, skill, map_summary)."),
      class = "dl_diagnostics_read_error"
    ))
  }
  out
}

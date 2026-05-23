# Shared helpers for apps/4a-distributional-learning-fit and
# apps/4b-distributional-learning-diagnostics.

`%||%` <- function(x, y) if (!is.null(x)) x else y

DL_YNAME <- "runoff_hourly"
DL_XNAMES <- c("rainfall_hourly", "snowmelt_hourly")

dl_default_cfg <- function() {
  list(
    yname = DL_YNAME,
    xnames = as.list(DL_XNAMES),
    na_action = "drop",
    min_obs = 5L,
    ntree = 500L,
    mtry = 1L,
    nodesize = 5L,
    nthreads = 1L
  )
}

dl_formula_cfg <- function() {
  list(
    yname = DL_YNAME,
    xnames = as.list(DL_XNAMES)
  )
}

dl_read_meta_full <- function(path) {
  if (!file.exists(path)) {
    return(list(dl_rqforest = dl_default_cfg()))
  }
  yaml::read_yaml(path)
}

dl_read_cfg <- function(path) {
  meta <- dl_read_meta_full(path)
  cfg <- meta$dl_rqforest
  if (is.null(cfg)) {
    return(dl_default_cfg())
  }
  cfg$xnames <- unlist(cfg$xnames, use.names = FALSE)
  cfg
}

dl_write_cfg <- function(path, cfg) {
  meta <- dl_read_meta_full(path)
  meta$dl_rqforest <- cfg
  yaml::write_yaml(meta, path)
}

dl_cfg_from_input <- function(input) {
  c(
    dl_formula_cfg(),
    list(
      na_action = input$na_action,
      min_obs = as.integer(max(1L, input$min_obs)),
      ntree = as.integer(max(10L, input$ntree)),
      mtry = as.integer(max(1L, min(2L, input$mtry))),
      nodesize = as.integer(max(1L, input$nodesize)),
      nthreads = as.integer(max(1L, input$nthreads))
    )
  )
}

dl_cfg_to_fit_args <- function(cfg) {
  c(
    list(
      yname = cfg$yname,
      xnames = unlist(cfg$xnames, use.names = FALSE),
      na_action = cfg$na_action %||% "drop",
      min_obs = as.integer(cfg$min_obs %||% 5L)
    ),
    cfg[setdiff(names(cfg), c("yname", "xnames", "na_action", "min_obs"))]
  )
}

dl_run_script4 <- function(root) {
  rscript <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
  script <- normalizePath(file.path(root, "scripts", "4-distributional_learning.r"), mustWork = FALSE)
  if (!file.exists(script)) {
    return(list(ok = FALSE, message = paste("Script not found:", script), output = character()))
  }
  owd <- getwd()
  on.exit(setwd(owd), add = TRUE)
  setwd(root)
  cmd <- paste(shQuote(normalizePath(rscript, mustWork = TRUE)), "--vanilla", shQuote(script))
  out <- tryCatch(
    suppressWarnings(system(cmd, intern = TRUE)),
    error = function(e) e
  )
  if (inherits(out, "error")) {
    return(list(ok = FALSE, message = conditionMessage(out), output = character()))
  }
  status <- attr(out, "status")
  ok <- is.null(status) || identical(as.integer(status), 0L)
  msg_out <- paste(as.character(out), collapse = "\n")
  list(
    ok = ok,
    message = if (!ok) paste("Script exited with status", status) else "Finished.",
    output = msg_out,
    status = status
  )
}

dl_param_hint <- function(text) {
  shiny::helpText(
    style = "margin-top: -10px; margin-bottom: 14px; font-size: 0.88em; color: #333; line-height: 1.4;",
    text
  )
}

dl_pal <- rev(c("#ff595e", "#ffca3a", "#8ac926", "#1982c4", "#6a4c93"))

dl_nearest_cell <- function(lon_click, lat_click, cells_tbl) {
  dx <- cells_tbl$y - lon_click
  dy <- cells_tbl$x - lat_click
  idx <- which.min(dx^2 + dy^2)
  cells_tbl$cell_id[idx]
}

dl_tile_dims <- function(xy_tbl, x_col = "y", y_col = "x") {
  ux <- sort(unique(xy_tbl[[x_col]]))
  uy <- sort(unique(xy_tbl[[y_col]]))
  w <- if (length(ux) > 1) stats::median(diff(ux)) else 0.25
  h <- if (length(uy) > 1) stats::median(diff(uy)) else 0.25
  c(width = w, height = h)
}

dl_compute_pp_long <- function(peak_hour_distributions) {
  dl_pp_long(peak_hour_distributions)
}

dl_compute_skill_scores <- function(peak_hour_distributions, dat) {
  dl_skill_scores(peak_hour_distributions, dat)
}

dl_load_diagnostics <- function(
  repo_root,
  peaks,
  dl_models = NULL,
  diagnostics_path = NULL,
  predictions_path = NULL
) {
  diagnostics_path <- diagnostics_path %||%
    fs::path(repo_root, "derived", dl_diagnostics_basename())
  predictions_path <- predictions_path %||%
    fs::path(repo_root, "derived", "era5_land_hourly_alps_dl_predictions.rds")

  saved <- dl_read_diagnostics(diagnostics_path)
  if (!inherits(saved, "dl_diagnostics_read_error")) {
    return(list(
      ok = TRUE,
      source = "saved",
      pp = saved$pp,
      skill = saved$skill,
      map_summary = saved$map_summary
    ))
  }

  dl_predictions <- dl_read_predictions_rds(predictions_path)
  rebuilt_predictions <- FALSE
  predictions_read_err <- NULL
  if (inherits(dl_predictions, "dl_predictions_read_error")) {
    predictions_read_err <- dl_predictions$message
    if (is.null(dl_models)) {
      return(list(
        ok = FALSE,
        source = NA_character_,
        message = paste(
          "Missing diagnostics RDS and predictions RDS;",
          "fitted models required to rebuild."
        )
      ))
    }
    dl_predictions <- dl_rebuild_predictions(peaks, dl_models)
    rebuilt_predictions <- TRUE
  }

  diag <- dl_build_diagnostics(dl_predictions, peaks)
  list(
    ok = TRUE,
    source = if (rebuilt_predictions) "rebuilt" else "recomputed",
    pp = diag$pp,
    skill = diag$skill,
    map_summary = diag$map_summary,
    predictions_error = predictions_read_err
  )
}

dl_fit_cell_preview <- function(peaks, cell_id, cfg) {
  cell_id <- as.integer(cell_id)
  cell_data <- peaks |> dplyr::filter(.data$cell_id == cell_id)
  if (nrow(cell_data) == 0L) {
    return(list(ok = FALSE, message = "No POT peaks for this cell."))
  }

  model <- do.call(dl_rqforest, c(list(data = cell_data), dl_cfg_to_fit_args(cfg)))
  x_df <- cell_data |> dplyr::select(dplyr::all_of(DL_XNAMES))
  dists <- stats::predict(model, newdata = x_df)

  preds <- cell_data |>
    dplyr::mutate(
      distribution_forest = dists,
      distribution_gp = purrr::map(.data$distribution_forest, convert_emp_to_gp)
    )

  pp_tbl <- dl_compute_pp_long(preds)
  skill_tbl <- dl_compute_skill_scores(preds, cell_data)
  skill_median <- stats::median(skill_tbl$skill_score, na.rm = TRUE)

  list(
    ok = TRUE,
    cell_id = cell_id,
    model = model,
    pp_tbl = pp_tbl,
    skill_tbl = skill_tbl,
    skill_median = skill_median,
    n_peaks = nrow(cell_data),
    fitted_at = Sys.time()
  )
}

dl_rebuild_predictions <- function(peaks, models_tbl) {
  peak_nests <- peaks |>
    tidyr::nest(.by = c("cell_id", "x", "y"))

  models_tbl |>
    dplyr::select(cell_id, x, y, dl_rqforest) |>
    dplyr::left_join(peak_nests, by = c("cell_id", "x", "y")) |>
    dplyr::mutate(
      distribution_forest = purrr::map2(
        .data$dl_rqforest,
        .data$data,
        \(m, d) stats::predict(m, newdata = d)
      )
    ) |>
    tidyr::unnest(c(.data$data, .data$distribution_forest)) |>
    dplyr::mutate(
      distribution_gp = purrr::map(.data$distribution_forest, convert_emp_to_gp)
    ) |>
    dplyr::select(
      cell_id,
      x,
      y,
      date,
      rainfall_hourly,
      snowmelt_hourly,
      runoff_hourly,
      dplyr::contains("distribution")
    )
}

dl_get_cell_model <- function(cell_id, dl_models) {
  if (is.null(dl_models)) {
    return(NULL)
  }
  cell_id <- as.integer(cell_id)
  row <- dl_models |> dplyr::filter(.data$cell_id == cell_id)
  if (nrow(row) != 1L) {
    return(NULL)
  }
  row$dl_rqforest[[1L]]
}

dl_read_predictions_rds <- function(path) {
  tryCatch(
    readRDS(path),
    error = function(e) {
      structure(
        list(message = conditionMessage(e)),
        class = "dl_predictions_read_error"
      )
    }
  )
}

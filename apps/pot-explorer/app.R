# POT selection explorer: hourly runoff, threshold, and peaks for one cell/year.
# Run from the package root:
#   shiny::runApp("apps/pot-explorer")
#
# Requires: hourly and POT RDS files from scripts 2 and 3 in data/; optional
#   inputs/pot_metadata.yaml (created with defaults by script 3 if missing).
# Suggests: yaml (for reading/writing pot_metadata.yaml)

library(shiny)
library(tidyverse)
library(sf)
library(fs)
library(rnaturalearth)
library(yaml)

repo_root <- here::here()

# Hourly timestamps are built as Jan 1 + hours; non-leap annual files spill a few
# rows into the next calendar year. That inflates max(year(date)) by 1. Here we tag
# each row with explorer_year: calendar year, except sparse tail years whose row
# count is tiny vs a typical year (detected below) — those rows attach to year-1.
explorer_year_from_dates <- function(tbl, date_col = "date") {
  dc <- tbl |>
    transmute(y = year(.data[[date_col]])) |>
    count(y, name = "n")
  med <- stats::median(dc$n)
  spill_years <- dc$y[dc$n < max(100L, as.integer(med * 0.01))]
  mutate(
    tbl,
    explorer_year = as.integer(
      if_else(
        year(.data[[date_col]]) %in% spill_years,
        year(.data[[date_col]]) - 1L,
        year(.data[[date_col]])
      )
    )
  )
}

join_peaks_explorer_year <- function(peaks_tbl, hourly_tbl) {
  lk <- hourly_tbl |> distinct(date, explorer_year)
  peaks_tbl |>
    dplyr::select(-dplyr::any_of("explorer_year")) |>
    left_join(lk, by = "date")
}

meta_path <- path(repo_root, "inputs", "pot_metadata.yaml")
hourly_path <- path(repo_root, "data", "era5_land_hourly_alps_all.rds")
peaks_path <- path(repo_root, "data", "era5_land_hourly_alps_peaks.rds")
threshold_path <- path(repo_root, "data", "era5_land_hourly_alps_pot_thresholds.rds")
read_pot_meta <- function(path) {
  if (!file.exists(path)) {
    return(list(quantile = 0.995, min_gap = 6L))
  }
  y <- read_yaml(path)
  q <- y$quantile
  mg <- y$min_gap
  list(
    quantile = if (is.null(q)) 0.995 else as.numeric(q),
    min_gap = if (is.null(mg)) 6L else as.integer(max(0, round(as.numeric(mg))))
  )
}

write_pot_meta <- function(path, quantile, min_gap) {
  write_yaml(
    list(quantile = as.numeric(quantile), min_gap = as.integer(min_gap)),
    path
  )
}

run_script3 <- function(root) {
  rscript <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
  script <- normalizePath(file.path(root, "scripts", "3-pot_spatial_eo.r"), mustWork = FALSE)
  if (!file.exists(script)) {
    return(list(
      ok = FALSE,
      message = paste("Script not found:", script)
    ))
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

nearest_cell <- function(lon_click, lat_click, cells_tbl) {
  dx <- cells_tbl$y - lon_click
  dy <- cells_tbl$x - lat_click
  idx <- which.min(dx^2 + dy^2)
  cells_tbl$cell_id[idx]
}

tile_dims <- function(xy_tbl, x_col = "y", y_col = "x") {
  ux <- sort(unique(xy_tbl[[x_col]]))
  uy <- sort(unique(xy_tbl[[y_col]]))
  w <- if (length(ux) > 1) stats::median(diff(ux)) else 0.25
  h <- if (length(uy) > 1) stats::median(diff(uy)) else 0.25
  c(width = w, height = h)
}

pot_meta0 <- read_pot_meta(meta_path)

hourly_all <- explorer_year_from_dates(read_rds(hourly_path))
peaks_all_init <- join_peaks_explorer_year(read_rds(peaks_path), hourly_all)
threshold_init <- read_rds(threshold_path)

cells_ref <- hourly_all |>
  distinct(cell_id, x, y) |>
  arrange(cell_id)

map_xlim <- range(cells_ref$y, na.rm = TRUE) + c(-0.5, 0.5)
map_ylim <- range(cells_ref$x, na.rm = TRUE) + c(-0.5, 0.5)

map_bbox <- st_bbox(
  c(xmin = map_xlim[1], xmax = map_xlim[2], ymin = map_ylim[1], ymax = map_ylim[2]),
  crs = st_crs(4326)
)

world_map <- ne_countries(scale = 50, returnclass = "sf") |>
  st_crop(map_bbox)

td <- tile_dims(cells_ref)
years_avail <- sort(unique(hourly_all$explorer_year))

ui <- fluidPage(
  titlePanel("Peaks over threshold (POT) explorer"),
  helpText(
    "Threshold comes from ",
    code("get_pot_events()"),
    " (quantile over ",
    code("runoff_hourly"),
    " per cell). ",
    "Adjust quantile / min gap, then ",
    strong("Re-run POT extraction"),
    " to rebuild peaks (runs ",
    code("scripts/3-pot_spatial_eo.r"),
    "). ",
    "Click the map to pick a cell."
  ),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      numericInput(
        "quantile",
        "Threshold quantile (per cell)",
        value = pot_meta0$quantile,
        min = 0,
        max = 1,
        step = 0.001
      ),
      numericInput(
        "min_gap",
        "Min gap between peaks (observations)",
        value = pot_meta0$min_gap,
        min = 0,
        step = 1
      ),
      actionButton("rerun_pot", "Re-run POT extraction", class = "btn-primary"),
      verbatimTextOutput("rerun_log"),
      hr(),
      selectInput("year", "Year", choices = years_avail, selected = max(years_avail)),
      selectInput(
        "cell_id",
        "Cell ID",
        choices = cells_ref$cell_id,
        selected = min(cells_ref$cell_id)
      ),
      verbatimTextOutput("cell_meta")
    ),
    mainPanel(
      width = 9,
      plotOutput("map_tiles", height = "360px", click = "map_click"),
      plotOutput("ts_runoff", height = "380px")
    )
  )
)

server <- function(input, output, session) {
  peaks_all <- reactiveVal(peaks_all_init)
  threshold_by_cell <- reactiveVal(threshold_init)
  rerun_log_txt <- reactiveVal("Adjust settings and click Re-run to rebuild POT outputs.")

  output$rerun_log <- renderText({
    rerun_log_txt()
  })

  observeEvent(input$rerun_pot, {
    q <- input$quantile
    mg <- input$min_gap
    if (!is.finite(q) || q <= 0 || q >= 1) {
      showNotification("Quantile must be strictly between 0 and 1.", type = "warning")
      return(invisible(NULL))
    }
    if (!is.finite(mg) || mg < 0) {
      showNotification("Min gap must be non-negative.", type = "warning")
      return(invisible(NULL))
    }

    write_pot_meta(meta_path, q, mg)

    withProgress(message = "Running POT script…", value = 0.5, {
      res <- run_script3(repo_root)
    })

    out_txt <- res$output
    if (!nzchar(out_txt)) {
      out_txt <- res$message
    } else {
      out_txt <- paste(out_txt, res$message, sep = "\n")
    }

    if (!isTRUE(res$ok)) {
      rerun_log_txt(paste0("ERROR\n", out_txt))
      showNotification(
        paste("POT script failed:", res$message),
        type = "error",
        duration = NULL
      )
      return(invisible(NULL))
    }

    peaks_all(join_peaks_explorer_year(read_rds(peaks_path), hourly_all))
    threshold_by_cell(read_rds(threshold_path))

    rerun_log_txt(paste(
      "Last run: OK",
      format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
      paste0("(wrote ", basename(meta_path), ")"),
      "",
      out_txt,
      sep = "\n"
    ))

    showNotification("POT extraction finished; data reloaded.", type = "message")
  })

  observeEvent(input$map_click, {
    cid <- nearest_cell(input$map_click$x, input$map_click$y, cells_ref)
    updateSelectInput(session, "cell_id", selected = cid)
  })

  output$cell_meta <- renderText({
    row <- cells_ref |> filter(cell_id == as.integer(input$cell_id))
    paste0(
      "lon = ", round(row$y, 3), ", lat = ", round(row$x, 3),
      "\n(runoff_hourly in mm)"
    )
  })

  thr_rv <- reactive({
    cid <- as.integer(input$cell_id)
    threshold_by_cell() |>
      filter(cell_id == cid) |>
      slice_head(n = 1) |>
      pull(threshold)
  })

  output$map_tiles <- renderPlot({
    ggplot(cells_ref, aes(y, x)) +
      geom_sf(
        data = world_map,
        inherit.aes = FALSE,
        fill = NA,
        linewidth = 1
      ) +
      geom_tile(
        aes(fill = cell_id == as.integer(input$cell_id)),
        colour = "grey30",
        linewidth = 0.3,
        width = td["width"],
        height = td["height"],
        alpha = 0.5
      ) +
      geom_text(aes(label = cell_id), size = 2.2, colour = "grey20") +
      scale_fill_manual(values = c(`TRUE` = "#ffd89b", `FALSE` = "#f7f7f7"), guide = "none") +
      coord_sf(xlim = map_xlim, ylim = map_ylim, expand = FALSE) +
      labs(x = "Longitude", y = "Latitude", title = "Cells (click to select)") +
      theme_minimal() +
      theme(panel.grid = element_blank())
  })

  output$ts_runoff <- renderPlot({
    yr <- as.integer(input$year)
    cid <- as.integer(input$cell_id)

    h <- hourly_all |>
      filter(cell_id == cid, explorer_year == yr)

    peaks_yr <- peaks_all() |>
      filter(cell_id == cid, explorer_year == yr)

    validate(
      need(nrow(h) > 0, "No hourly data for this cell and year — pick another combination.")
    )

    thr <- thr_rv()

    ggplot(h, aes(date, runoff_hourly)) +
      geom_line(linewidth = 0.2, colour = "#2c7fb8") +
      geom_hline(
        yintercept = thr,
        colour = "#d95f0e",
        linetype = "dashed",
        linewidth = 0.8
      ) +
      geom_point(
        data = peaks_yr,
        mapping = aes(date, runoff_hourly),
        shape = 21,
        size = 3,
        stroke = 1,
        fill = "#d95f0e",
        colour = "white"
      ) +
      theme_bw() +
      labs(
        title = paste0("Hourly runoff | cell ", cid, ", year ", yr),
        x = NULL,
        y = "Hourly runoff (mm)"
      )
  })
}

shinyApp(ui, server)

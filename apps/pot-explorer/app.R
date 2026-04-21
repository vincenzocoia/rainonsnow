# POT selection explorer: hourly runoff, threshold, and peaks for one cell/year.
# Run from the package root:
#   shiny::runApp("apps/pot-explorer")
#
# Requires: hourly and POT RDS files from scripts 2 and 3 in data/.

library(shiny)
library(tidyverse)
library(sf)
library(fs)
library(rnaturalearth)

repo_root <- here::here()
hourly_path <- path(repo_root, "data", "era5_land_hourly_alps_all.rds")
peaks_path <- path(repo_root, "data", "era5_land_hourly_alps_peaks.rds")
threshold_path <- path(repo_root, "data", "era5_land_hourly_alps_pot_thresholds.rds")

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

hourly_all <- read_rds(hourly_path)
peaks_all <- read_rds(peaks_path)
threshold_by_cell <- read_rds(threshold_path)

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
years_avail <- sort(unique(year(hourly_all$date)))

ui <- fluidPage(
  titlePanel("Peaks over threshold (POT) explorer"),
  helpText(
    "Threshold comes from ",
    code("get_pot_events()"),
    " (column ",
    code("threshold"),
    " in the POT file when present). ",
    "Click the map to pick a cell."
  ),
  sidebarLayout(
    sidebarPanel(
      width = 3,
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
    threshold_by_cell |>
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
      filter(cell_id == cid, year(date) == yr)

    peaks_yr <- peaks_all |>
      filter(cell_id == cid, year(date) == yr)

    validate(
      need(nrow(h) > 0, "No hourly data for this cell and year — pick another combination.")
    )

    thr <- thr_rv()
    thr_lbl <- if (is.finite(thr)) sprintf("%.3g mm", thr) else "NA"

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

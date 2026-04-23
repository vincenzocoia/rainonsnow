# Runoff marginal explorer: DL mixture return curves vs naive POT-only marginals.
# Run from the package root:
#   shiny::runApp("apps/runoff-marginals-explorer")
#
# Requires:
#   - data/era5_land_hourly_alps_peaks.rds
#   - data/era5_land_hourly_alps_dl_return_levels.rds (scripts/5-runoff_marginals.r)
#
# Suggests: shiny, tidyverse, sf, rnaturalearth, distionary, famish

library(shiny)
library(tidyverse)
library(sf)
library(fs)
library(rnaturalearth)
library(distionary)
library(famish)

repo_root <- here::here()
devtools::load_all(repo_root, quiet = TRUE)

peaks_path <- path(repo_root, "data", "era5_land_hourly_alps_peaks.rds")
marginal_levels_path <- path(repo_root, "data", "era5_land_hourly_alps_dl_return_levels.rds")

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

read_marginal_return_levels <- function(path) {
  if (!file.exists(path)) {
    return(NULL)
  }
  raw <- readRDS(path)
  if (!is.data.frame(raw)) {
    return(NULL)
  }
  if ("return_period" %in% names(raw) && !"return_period_years" %in% names(raw)) {
    raw <- dplyr::rename(raw, return_period_years = return_period)
  }
  if (!"model" %in% names(raw)) {
    return(NULL)
  }
  ## Script 5 used names_prefix = "return_" once, so model could stay
  ## "levels_forest" / "levels_gp" — those never matched scale_colour_manual.
  raw |>
    dplyr::mutate(
      model = dplyr::case_when(
        .data$model %in% c("forest", "levels_forest") ~ "Random Forest",
        .data$model %in% c("gp", "levels_gp") ~ "GP conversion",
        .data$model %in% c("Random Forest", "GP conversion") ~ .data$model,
        TRUE ~ .data$model
      )
    )
}

## Apply shared options to FM ggplot: optional log y, optional matched y limits.
apply_fm_y_scales <- function(plot, log_y, y_limits) {
  if (isTRUE(log_y)) {
    plot <- plot + ggplot2::scale_y_log10()
  }
  if (!is.null(y_limits) && length(y_limits) == 2L && all(is.finite(y_limits))) {
    if (isTRUE(log_y)) {
      y_limits <- pmax(y_limits, .Machine$double.eps)
    }
    plot <- plot + ggplot2::coord_cartesian(ylim = y_limits)
  }
  plot
}

# --- load static inputs ------------------------------------------------------

peaks_ok <- file.exists(peaks_path)
dat <- if (peaks_ok) readRDS(peaks_path) else NULL

marginal_long <- read_marginal_return_levels(marginal_levels_path)
marginal_ok <- is.data.frame(marginal_long) && nrow(marginal_long) > 0L

cells_ref <- if (!is.null(dat)) {
  dat |>
    distinct(cell_id, x, y) |>
    arrange(cell_id)
} else {
  tibble(cell_id = integer(), x = numeric(), y = numeric())
}

map_xlim <- if (nrow(cells_ref) > 0) range(cells_ref$y, na.rm = TRUE) + c(-0.5, 0.5) else c(-1, 1)
map_ylim <- if (nrow(cells_ref) > 0) range(cells_ref$x, na.rm = TRUE) + c(-0.5, 0.5) else c(-1, 1)
map_bbox <- st_bbox(
  c(xmin = map_xlim[1], xmax = map_xlim[2], ymin = map_ylim[1], ymax = map_ylim[2]),
  crs = st_crs(4326)
)
world_map <- tryCatch(
  rnaturalearth::ne_countries(scale = 50, returnclass = "sf") |>
    st_crop(map_bbox),
  error = function(e) NULL
)
td <- if (nrow(cells_ref) > 0) tile_dims(cells_ref) else c(width = 0.25, height = 0.25)

rp_grid <- rp_reporting()

## Orange-red vs blue (Wong palette): easier to distinguish than green vs blue.
colour_forest <- "#D55E00"
colour_gp <- "#0072B2"
dl_colours <- c(
  `Random Forest` = colour_forest,
  `GP conversion` = colour_gp
)
naive_colours <- c(
  Empirical = colour_forest,
  `GP (naive)` = colour_gp
)

## Shared plot height so DL vs naive frequency–magnitude panels match exactly.
fm_plot_height <- "380px"
map_plot_height <- "260px"

ui <- fluidPage(
  titlePanel("Runoff marginals explorer"),
  uiOutput("data_banner"),
  helpText(
    "Below the map, two frequency–magnitude panels are shown side by side at the same size: ",
    "left — distributional learning (Random Forest vs GP conversion); ",
    "right — naive POT marginals (",
    code("dst_empirical"),
    " vs ",
    code("fit_dst_gp"),
    "). Same POT event axis; sidebar toggles sync y-axis limits across panels and apply log ",
    code("y"),
    " on both. Click the map or choose a cell."
  ),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      selectInput(
        "cell_id",
        "Cell ID",
        choices = if (nrow(cells_ref) > 0) cells_ref$cell_id else integer(),
        selected = if (nrow(cells_ref) > 0) min(cells_ref$cell_id) else NULL
      ),
      checkboxInput(
        "match_y_scales",
        "Match y-axis limits (DL vs naive)",
        value = FALSE
      ),
      checkboxInput(
        "log_y_return",
        "Log scale for return level (y), both panels",
        value = FALSE
      ),
      verbatimTextOutput("cell_meta")
    ),
    mainPanel(
      width = 9,
      fluidRow(
        column(12, plotOutput("map_tiles", height = map_plot_height, click = "map_click"))
      ),
      fluidRow(
        column(6, plotOutput("fm_dl", height = fm_plot_height)),
        column(6, plotOutput("fm_naive", height = fm_plot_height))
      )
    )
  )
)

server <- function(input, output, session) {
  observeEvent(input$map_click, {
    req(nrow(cells_ref) > 0)
    cid <- nearest_cell(input$map_click$x, input$map_click$y, cells_ref)
    updateSelectInput(session, "cell_id", selected = cid)
  })

  output$data_banner <- renderUI({
    msg <- character()
    if (!peaks_ok) {
      msg <- c(msg, basename(peaks_path))
    }
    if (!marginal_ok) {
      msg <- c(msg, basename(marginal_levels_path))
    }
    if (length(msg) == 0L) {
      return(NULL)
    }
    div(
      class = "alert alert-warning",
      style = "padding:8px;margin-bottom:12px;",
      strong("Missing data: "),
      paste(msg, collapse = ", "),
      " — run ",
      code("scripts/3-pot_spatial_eo.r"),
      ", ",
      code("scripts/4-distributional_learning.r"),
      ", ",
      code("scripts/5-runoff_marginals.r"),
      "."
    )
  })

  output$cell_meta <- renderText({
    req(peaks_ok, nrow(cells_ref) > 0)
    row <- cells_ref |> filter(cell_id == as.integer(input$cell_id))
    sub <- dat |>
      filter(cell_id == as.integer(input$cell_id)) |>
      filter(is.finite(runoff_hourly))
    nep <- if (nrow(sub) > 0) {
      yr <- as.integer(format(as.Date(sub$date), "%Y"))
      span <- max(yr, na.rm = TRUE) - min(yr, na.rm = TRUE) + 1L
      nrow(sub) / max(span, 1L)
    } else {
      NA_real_
    }
    paste0(
      "lon = ", round(row$y, 3), ", lat = ", round(row$x, 3),
      "\nPOT peaks: ", nrow(sub),
      if (is.finite(nep)) paste0("; ~", signif(nep, 4), " peaks/year") else "",
      "\n(runoff_hourly in mm)"
    )
  })

  output$map_tiles <- renderPlot({
    req(!is.null(world_map), nrow(cells_ref) > 0)
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
      scale_fill_manual(values = c(`TRUE` = "#b3cde3", `FALSE` = "#f7f7f7"), guide = "none") +
      coord_sf(xlim = map_xlim, ylim = map_ylim, expand = FALSE) +
      labs(x = "Longitude", y = "Latitude", title = "Cells (click to select)") +
      theme_minimal() +
      theme(panel.grid = element_blank())
  })

  dl_cell_tbl <- reactive({
    if (!marginal_ok || is.null(marginal_long)) {
      return(tibble::tibble())
    }
    cid <- as.integer(input$cell_id)
    marginal_long |>
      filter(cell_id == cid) |>
      filter(is.finite(return_level))
  })

  peaks_cell <- reactive({
    req(peaks_ok, input$cell_id)
    dat |>
      filter(cell_id == as.integer(input$cell_id)) |>
      pull(runoff_hourly)
  })

  cell_nep <- reactive({
    req(peaks_ok, input$cell_id)
    sub <- dat |>
      filter(cell_id == as.integer(input$cell_id)) |>
      filter(is.finite(runoff_hourly))
    req(nrow(sub) > 0)
    yr <- as.integer(format(as.Date(sub$date), "%Y"))
    span <- max(yr, na.rm = TRUE) - min(yr, na.rm = TRUE) + 1L
    nrow(sub) / max(span, 1L)
  })

  naive_cell_tbl <- reactive({
    if (!peaks_ok) {
      return(NULL)
    }
    cid <- as.integer(input$cell_id)
    sub <- dat |>
      filter(cell_id == cid) |>
      filter(is.finite(runoff_hourly))
    if (nrow(sub) < 2L) {
      return(NULL)
    }
    y <- sub$runoff_hourly[is.finite(sub$runoff_hourly)]
    if (length(y) < 2L) {
      return(NULL)
    }
    yr <- as.integer(format(as.Date(sub$date), "%Y"))
    span <- max(yr, na.rm = TRUE) - min(yr, na.rm = TRUE) + 1L
    nep <- nrow(sub) / max(span, 1L)

    emp <- tryCatch(distionary::dst_empirical(y), error = function(e) NULL)
    gp <- tryCatch(famish::fit_dst_gp(y), error = function(e) NULL)
    if (is.null(emp) || is.null(gp)) {
      return(NULL)
    }

    fm_emp <- tryCatch(
      enframe_at_events(emp, nep, return_periods = rp_grid),
      error = function(e) NULL
    )
    fm_gp <- tryCatch(
      enframe_at_events(gp, nep, return_periods = rp_grid),
      error = function(e) NULL
    )
    if (is.null(fm_emp) || is.null(fm_gp)) {
      return(NULL)
    }

    dplyr::bind_rows(
      dplyr::mutate(fm_emp, model = "Empirical"),
      dplyr::mutate(fm_gp, model = "GP (naive)")
    ) |>
      dplyr::rename(return_level = `return`) |>
      dplyr::filter(is.finite(return_level))
  })

  shared_y_limits <- reactive({
    if (!isTRUE(input$match_y_scales)) {
      return(NULL)
    }
    dl <- dl_cell_tbl()
    nv <- naive_cell_tbl()
    if (nrow(dl) == 0L || is.null(nv) || nrow(nv) == 0L) {
      return(NULL)
    }
    vals <- c(dl$return_level, nv$return_level)
    vals <- vals[is.finite(vals)]
    log_y <- isTRUE(input$log_y_return)
    if (log_y) {
      vals <- vals[vals > 0]
    }
    if (length(vals) == 0L) {
      return(NULL)
    }
    rng <- range(vals, na.rm = TRUE)
    if (log_y) {
      c(max(rng[1] / 1.08, .Machine$double.eps), rng[2] * 1.08)
    } else {
      pad <- diff(rng) * 0.06
      if (!is.finite(pad) || pad <= 0) {
        pad <- max(abs(rng[2]), abs(rng[1]), 1e-6) * 0.06
      }
      c(rng[1] - pad, rng[2] + pad)
    }
  })

  output$fm_dl <- renderPlot({
    validate <- shiny::validate
    need <- shiny::need
    validate(need(marginal_ok, "Load marginal return levels (script 5)."))
    cid <- as.integer(input$cell_id)
    sub <- dl_cell_tbl()
    validate(need(nrow(sub) > 0, "No DL marginal rows for this cell."))

    log_y <- isTRUE(input$log_y_return)
    y_lab <- if (log_y) "Return level (mm, log scale)" else "Return level (mm)"
    ylim_use <- if (isTRUE(input$match_y_scales)) shared_y_limits() else NULL

    p <- ggplot(sub, aes(return_period_years, return_level, colour = model, group = model)) +
      geom_line(linewidth = 0.95, alpha = 0.9) +
      geom_point(size = 1.8, alpha = 0.85) +
      scale_x_log10(breaks = rp_grid, minor_breaks = NULL) +
      scale_colour_manual(
        values = dl_colours,
        breaks = names(dl_colours),
        limits = names(dl_colours)
      ) +
      theme_bw() +
      labs(
        title = paste0("DL marginals | cell ", cid),
        subtitle = "Mixture of hourly predictive distributions (POT event axis)",
        x = "Return period (years)",
        y = y_lab,
        colour = NULL
      ) +
      theme(legend.position = "bottom")

    apply_fm_y_scales(p, log_y, ylim_use)
  })

  output$fm_naive <- renderPlot({
    validate <- shiny::validate
    need <- shiny::need
    validate(need(peaks_ok, "Peaks RDS missing."))
    combined <- naive_cell_tbl()
    validate(
      need(!is.null(combined) && nrow(combined) > 0, paste(
        "Naive curves unavailable:",
        "need ≥2 POT peaks and successful empirical + GP fits."
      ))
    )

    cid <- as.integer(input$cell_id)
    nep <- cell_nep()
    log_y <- isTRUE(input$log_y_return)
    y_lab <- if (log_y) "Return level (mm, log scale)" else "Return level (mm)"
    ylim_use <- if (isTRUE(input$match_y_scales)) shared_y_limits() else NULL

    p <- ggplot(combined, aes(return_period_years, return_level, colour = model, group = model)) +
      geom_line(linewidth = 0.95, alpha = 0.9) +
      geom_point(size = 1.8, alpha = 0.85) +
      scale_x_log10(breaks = rp_grid, minor_breaks = NULL) +
      scale_colour_manual(
        values = naive_colours,
        breaks = names(naive_colours),
        limits = names(naive_colours)
      ) +
      theme_bw() +
      labs(
        title = paste0("Naive POT marginals | cell ", cid),
        subtitle = paste0(
          "dst_empirical vs fit_dst_gp on raw peaks; ~",
          signif(nep, 4),
          " peaks/year"
        ),
        x = "Return period (years)",
        y = y_lab,
        colour = NULL
      ) +
      theme(legend.position = "bottom")

    apply_fm_y_scales(p, log_y, ylim_use)
  })
}

shinyApp(ui, server)

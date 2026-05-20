# Return-period explorer: map of marginal return levels and cell diagnostics.
# Run from the package root:
#   shiny::runApp("apps/return-level-explorer")
#
# Requires:
#   - derived/era5_land_hourly_alps_peaks.rds
#   - derived/era5_land_hourly_alps_dl_rqforest_models.rds  (script 4)
#   - derived/era5_land_hourly_alps_dl_return_levels.rds  (script 5 flat table; preferred)
#     OR derived/era5_land_hourly_alps_dl_marginal_return_levels.rds  (bundle with marginal_return_levels_long)
#     OR derived/era5_land_hourly_alps_dl_predictions.rds as slow fallback (mixtures built in app)

library(shiny)
library(tidyverse)
library(sf)
library(fs)
library(rnaturalearth)
library(probaverse)

repo_root <- here::here()
devtools::load_all(repo_root, quiet = TRUE)

peaks_path <- path(repo_root, "derived", "era5_land_hourly_alps_peaks.rds")
dl_path <- path(repo_root, "derived", "era5_land_hourly_alps_dl_predictions.rds")
models_path <- path(repo_root, "derived", "era5_land_hourly_alps_dl_rqforest_models.rds")
marginal_flat_path <- path(repo_root, "derived", "era5_land_hourly_alps_dl_return_levels.rds")
marginal_bundle_path <- path(repo_root, "derived", "era5_land_hourly_alps_dl_marginal_return_levels.rds")

## Script 5 writes a long table with model labels "Random Forest" / "GP conversion";
## this app plots using the names produced by its slow-path rebuild.
normalize_marginal_levels_for_app <- function(tbl) {
  if (!is.data.frame(tbl) || nrow(tbl) == 0L) {
    return(tibble::tibble())
  }
  out <- tbl
  if ("return_period" %in% names(out) && !"return_period_years" %in% names(out)) {
    out <- dplyr::rename(out, return_period_years = return_period)
  }
  out |>
    dplyr::mutate(
      model = dplyr::case_when(
        .data$model %in% c("Random Forest", "forest", "Forest mixture (empirical)") ~
          "Forest mixture (empirical)",
        .data$model %in% c("GP conversion", "gp", "Forest mixture + GP tail") ~
          "Forest mixture + GP tail",
        TRUE ~ .data$model
      )
    )
}

## Rain–snow panel only: grid + forecasts for one cell (was all cells at startup — very slow).
cell_rain_snow_prep_one <- function(cid, dat_all, dl_models_tbl) {
  cid <- as.integer(cid)
  sub <- dat_all |> dplyr::filter(.data$cell_id == cid)
  if (nrow(sub) == 0L) {
    return(NULL)
  }
  row_dl <- dl_models_tbl |> dplyr::filter(.data$cell_id == cid)
  if (nrow(row_dl) != 1L) {
    return(NULL)
  }
  g <- grid_from_scatter(
    sub$rainfall_hourly,
    sub$snowmelt_hourly,
    data = sub,
    mult = 1.3
  )
  g <- dplyr::rename(g, rainfall_hourly = x, snowmelt_hourly = y)
  forecast_forest <- stats::predict(row_dl$dl_rqforest[[1]], newdata = g)
  forecast_gp <- purrr::map(
    forecast_forest,
    \(d) tryCatch(convert_emp_to_gp(d), error = function(e) d)
  )
  list(data = sub, grid = g, forecast_gp = forecast_gp)
}

pal <- rev(c("#ff595e", "#ffca3a", "#8ac926", "#1982c4", "#6a4c93"))

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

dat <- read_rds(peaks_path)

cells_ref <- dat |>
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

return_periods <- c(2, 5, 10, 20, 50, 100, 200, 500, 1000)

data_ready <- file.exists(models_path) &&
  (
    file.exists(marginal_flat_path) ||
      file.exists(marginal_bundle_path) ||
      file.exists(dl_path)
  )

ui <- fluidPage(
  titlePanel("Marginal return levels and rain–snowmelt drivers"),
  uiOutput("data_banner"),
  helpText(
    "Select a return period to colour the map by marginal hourly runoff return level ",
    "(mixture of predictive distributions per peak hour; GP-adjusted tail). ",
    "Click a cell for frequency–magnitude curves (forest mixture vs GP tail) and a rain–snowmelt ",
    "likelihood surface for the runoff depth at the chosen return period ",
    "(event frequency accounts for POT events per year)."
  ),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      selectInput(
        "return_period",
        "Return period (years)",
        choices = return_periods,
        selected = 200
      ),
      selectInput(
        "cell_id",
        "Cell ID",
        choices = cells_ref$cell_id,
        selected = min(cells_ref$cell_id)
      ),
      verbatimTextOutput("cell_meta"),
      verbatimTextOutput("prep_status")
    ),
    mainPanel(
      width = 9,
      plotOutput("map_return", height = "340px", click = "map_click"),
      fluidRow(
        column(6, plotOutput("fm_curve", height = "320px")),
        column(6, plotOutput("rain_snow", height = "320px"))
      )
    )
  )
)

server <- function(input, output, session) {
  values <- reactiveValues(
    ready = FALSE,
    marginal_return_levels_long = NULL,
    dl_models_tbl = NULL,
    error = NULL
  )

  output$data_banner <- renderUI({
    if (data_ready) {
      return(NULL)
    }
    missing <- c()
    if (!file.exists(models_path)) {
      missing <- c(missing, basename(models_path))
    }
    if (
      !file.exists(marginal_flat_path) &&
        !file.exists(marginal_bundle_path) &&
        !file.exists(dl_path)
    ) {
      missing <- c(
        missing,
        paste0(
          basename(marginal_flat_path), ", ",
          basename(marginal_bundle_path), ", or ",
          basename(dl_path)
        )
      )
    }
    helpText(
      strong("Missing data: "),
      paste(missing, collapse = "; "),
      ". Run ",
      code("scripts/4-distributional_learning.r"),
      ". For faster startup, also run ",
      code("scripts/5-runoff_marginals.r"),
      " (precomputed marginal return levels). Then reload this app."
    )
  })

  observe({
    if (!data_ready) {
      return(invisible(NULL))
    }

    values$error <- NULL
    values$ready <- FALSE

    tryCatch(
      {
        withProgress(
          message = "Preparing maps and rain–snowmelt surfaces…",
          detail = "Reading DL outputs",
          value = 0,
          {
            dl_models <- read_rds(models_path)

            if (file.exists(marginal_bundle_path)) {
              incProgress(0.15, detail = "Reading precomputed marginal return levels (bundle)")
              bundle <- read_rds(marginal_bundle_path)
              marginal_return_levels_long <- bundle$marginal_return_levels_long
              rp_saved <- bundle$return_periods
              if (!identical(unname(rp_saved), unname(return_periods))) {
                showNotification(
                  paste(
                    "Return periods in the precomputed RDS differ from `return_periods` in app.R.",
                    "Align them or update the drop-down choices."
                  ),
                  type = "warning",
                  duration = 12
                )
              }
            } else if (file.exists(marginal_flat_path)) {
              incProgress(0.15, detail = "Reading script 5 marginal return levels")
              marginal_return_levels_long <- normalize_marginal_levels_for_app(
                read_rds(marginal_flat_path)
              )
              rp_saved <- sort(unique(marginal_return_levels_long$return_period_years))
              if (!identical(unname(rp_saved), unname(return_periods))) {
                showNotification(
                  paste(
                    "Return periods in ", basename(marginal_flat_path),
                    " differ from `return_periods` in app.R.",
                    "Align them or update the drop-down choices."
                  ),
                  type = "warning",
                  duration = 12
                )
              }
            } else {
              incProgress(0.05, detail = "Reading hourly predictions (slow path: marginal mixtures)")
              peak_hour_distributions <- read_rds(dl_path)

              num_pot_events <- dat |>
                group_by(cell_id, x, y) |>
                summarise(
                  num_events_per_year = n() / (diff(range(year(date))) + 1),
                  .groups = "drop"
                )

              marginals <- peak_hour_distributions |>
                group_by(cell_id, x, y) |>
                summarise(
                  marginal_runoff = list(mix2(distribution_forest, na_action_dst = "drop")),
                  .groups = "drop"
                )

              marginalsgp <- peak_hour_distributions |>
                group_by(cell_id, x, y) |>
                summarise(
                  marginal_runoff = list(mix2(distribution_gp, na_action_dst = "drop")),
                  .groups = "drop"
                )

              marginals <- marginals |>
                left_join(num_pot_events, by = c("cell_id", "x", "y"))
              marginalsgp <- marginalsgp |>
                left_join(num_pot_events, by = c("cell_id", "x", "y"))

              forest_return_levels <- marginals |>
                mutate(
                  df = map2(
                    marginal_runoff,
                    num_events_per_year,
                    enframe_at_events,
                    return_periods = return_periods
                  )
                ) |>
                select(cell_id, x, y, df) |>
                unnest(df)

              gp_return_levels <- marginalsgp |>
                mutate(
                  df = map2(
                    marginal_runoff,
                    num_events_per_year,
                    enframe_at_events,
                    return_periods = return_periods
                  )
                ) |>
                select(cell_id, x, y, df) |>
                unnest(df)

              return_levels_joined <- left_join(
                forest_return_levels,
                gp_return_levels |> select(cell_id, x, y, return_period, return),
                by = c("cell_id", "x", "y", "return_period"),
                suffix = c("_emp", "_gp")
              )

              marginal_return_levels_long <- return_levels_joined |>
                pivot_longer(
                  c(return_emp, return_gp),
                  names_to = "model",
                  values_to = "return_level",
                  names_prefix = "return_"
                ) |>
                mutate(
                  model = case_when(
                    model == "emp" ~ "Forest mixture (empirical)",
                    model == "gp" ~ "Forest mixture + GP tail",
                    TRUE ~ model
                  )
                )
              incProgress(0.55, detail = "Marginal mixtures finished")
            }

            incProgress(0.95, detail = "Done")
            values$marginal_return_levels_long <- marginal_return_levels_long
            values$dl_models_tbl <- dl_models
            values$ready <- TRUE
          }
        )
      },
      error = function(e) {
        values$error <- conditionMessage(e)
        values$ready <- FALSE
        showNotification(paste("Prep failed:", values$error), type = "error", duration = NULL)
      }
    )
  })

  output$prep_status <- renderText({
    if (!data_ready) {
      return("Waiting for script 4 outputs.")
    }
    err <- values$error
    if (length(err) && nzchar(err)) {
      return(paste("Error:", err))
    }
    if (!values$ready) {
      return("Loading…")
    }
    "Ready."
  })

  ## Deferred: rain–snow likelihood surface is built only for the selected cell.
  rain_snow_prep <- reactive({
    req(values$ready, values$dl_models_tbl)
    cid <- as.integer(input$cell_id)
    cell_rain_snow_prep_one(cid, dat, values$dl_models_tbl)
  })

  map_df <- reactive({
    req(values$ready, values$marginal_return_levels_long)
    rp <- as.numeric(input$return_period)
    values$marginal_return_levels_long |>
      filter(return_period_years == rp, model == "Forest mixture + GP tail") |>
      select(cell_id, x, y, return_level)
  })

  observeEvent(input$map_click, {
    cid <- nearest_cell(input$map_click$x, input$map_click$y, cells_ref)
    updateSelectInput(session, "cell_id", selected = cid)
  })

  output$cell_meta <- renderText({
    row <- cells_ref |> filter(cell_id == as.integer(input$cell_id))
    paste0(
      "lon = ", round(row$y, 3), ", lat = ", round(row$x, 3),
      "\nPOT peaks (hourly runoff, mm)"
    )
  })

  output$map_return <- renderPlot({
    req(values$ready)
    rp <- as.numeric(input$return_period)
    d <- map_df() |>
      mutate(selected = cell_id == as.integer(input$cell_id))

    ggplot(d, aes(y, x)) +
      geom_sf(
        data = world_map,
        inherit.aes = FALSE,
        fill = NA,
        linewidth = 1
      ) +
      geom_tile(
        aes(fill = return_level),
        linewidth = 0.25,
        colour = scales::alpha("grey35", 0.35),
        width = td["width"],
        height = td["height"],
        alpha = 0.88
      ) +
      geom_tile(
        data = dplyr::filter(d, selected),
        aes(y, x),
        inherit.aes = FALSE,
        fill = NA,
        colour = "grey10",
        linewidth = 0.65,
        width = td["width"],
        height = td["height"]
      ) +
      geom_text(aes(label = cell_id), size = 2.2, colour = "grey15") +
      scale_fill_gradientn(
        paste0(rp, "-year return\nlevel (mm)"),
        colours = pal,
        na.value = "grey92"
      ) +
      coord_sf(xlim = map_xlim, ylim = map_ylim, expand = FALSE) +
      labs(
        x = "Longitude",
        y = "Latitude",
        title = paste0("Marginal return level (GP tail) at ", rp, "-year return period")
      ) +
      theme_minimal() +
      theme(panel.grid = element_blank())
  })

  output$fm_curve <- renderPlot({
    req(values$ready, values$marginal_return_levels_long)
    cid <- as.integer(input$cell_id)
    values$marginal_return_levels_long |>
      filter(cell_id == cid) |>
      ggplot(aes(return_period_years, return_level, colour = model)) +
      geom_line(linewidth = 0.9, alpha = 0.85) +
      geom_vline(
        xintercept = as.numeric(input$return_period),
        linetype = "dotted",
        colour = "grey35"
      ) +
      scale_x_log10("Return period (years)") +
      scale_colour_manual(values = c("black", "red3")) +
      labs(
        title = paste0("Frequency–magnitude | cell ", cid),
        y = "Return level (mm)",
        colour = NULL
      ) +
      theme_bw() +
      theme(legend.position = "bottom")
  })

  output$rain_snow <- renderPlot({
    req(values$ready)
    cid <- as.integer(input$cell_id)
    rp <- as.numeric(input$return_period)

    prep <- rain_snow_prep()
    validate(
      need(
        !is.null(prep),
        "No POT data or DL model row for this cell."
      )
    )

    mag_tbl <- values$marginal_return_levels_long |>
      filter(
        cell_id == cid,
        return_period_years == rp,
        model == "Forest mixture + GP tail"
      )
    validate(need(nrow(mag_tbl) == 1L, "Return level not available for this cell."))
    magnitude <- mag_tbl$return_level[1]

    g <- prep$grid
    fgp <- prep$forecast_gp
    validate(need(length(fgp) == nrow(g), "Forecast grid mismatch."))

    g <- g |>
      mutate(prob_responsible = distribution_likeliness(fgp, magnitude))

    pts <- prep$data |>
      select(rainfall_hourly, snowmelt_hourly)

    ggplot(g, aes(rainfall_hourly, snowmelt_hourly)) +
      geom_contour_filled(aes(z = prob_responsible), alpha = 0.75) +
      geom_contour(aes(z = prob_responsible), colour = scales::alpha("grey15", 0.35), linewidth = 0.2) +
      geom_point(
        data = pts,
        inherit.aes = FALSE,
        aes(rainfall_hourly, snowmelt_hourly),
        colour = "black",
        alpha = 0.12,
        size = 0.35
      ) +
      coord_cartesian(expand = FALSE) +
      scale_fill_viridis_c("Relative\ndensity", option = "C", end = 0.95) +
      labs(
        title = paste0(
          "Rain vs snowmelt | ", rp,
          "-year runoff (", signif(magnitude, 4), " mm)"
        ),
        x = "Rainfall (mm/h)",
        y = "Snowmelt (mm/h)"
      ) +
      theme_bw() +
      theme(legend.position = "right")
  })
}

shinyApp(ui, server)

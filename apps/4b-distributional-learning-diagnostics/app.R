# Distributional learning — diagnostics (script 4 outputs).
# Run from the package root:
#   shiny::runApp("apps/4b-distributional-learning-diagnostics")
#
# Requires: derived/era5_land_hourly_alps_peaks.rds and
#   era5_land_hourly_alps_dl_rqforest_models.rds (script 4).
# P–P / skill: era5_land_hourly_alps_dl_diagnostics.rds (script 4); falls back to
#   predictions RDS or rebuild from models if diagnostics are missing.

library(shiny)
library(tidyverse)
library(sf)
library(fs)
library(rnaturalearth)
library(probaverse)
library(distionary)

repo_root <- here::here()
devtools::load_all(repo_root, quiet = TRUE)
source(file.path(repo_root, "apps", "ros_theme.R"))
source(fs::path(repo_root, "apps", "dl_shared.R"))

peaks_path <- fs::path(repo_root, "derived", "era5_land_hourly_alps_peaks.rds")
dl_path <- fs::path(repo_root, "derived", "era5_land_hourly_alps_dl_predictions.rds")
diag_path <- fs::path(repo_root, "derived", dl_diagnostics_basename())
models_path <- fs::path(repo_root, "derived", "era5_land_hourly_alps_dl_rqforest_models.rds")

peaks_ok <- file.exists(peaks_path)
models_ok <- file.exists(models_path)
diagnostics_ok <- file.exists(diag_path)

ui <- fluidPage(
  theme = ros_bs_theme(),
  ros_header(
    "Distributional learning — diagnostics",
    "Calibration, skill against the marginal, and conditional runoff CDFs per cell."
  ),
  uiOutput("data_banner"),
  helpText(
    "Explore fitted models from script 4. Models load at startup (for CDF clicks on any cell). ",
    "Load P–P / skill / map from the diagnostics RDS (fast); otherwise recompute from predictions or models."
  ),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      verbatimTextOutput("load_status"),
      actionButton(
        "load_diagnostics",
        "Load P–P and skill diagnostics",
        class = "btn-primary",
        width = "100%"
      ),
      helpText(
        style = "font-size: 0.85em;",
        "Reads diagnostics RDS when present; otherwise recomputes (~1–2 min)."
      ),
      actionButton("reload_models", "Reload fitted models", class = "btn-default", width = "100%"),
      hr(),
      selectInput("cell_id", "Cell ID", choices = integer(0)),
      verbatimTextOutput("cell_meta")
    ),
    mainPanel(
      width = 9,
      plotOutput("map_skill", height = "320px", click = "map_click"),
      fluidRow(
        column(6, plotOutput("pp_plot", height = "300px")),
        column(6, plotOutput("skill_plot", height = "300px"))
      ),
      fluidRow(
        column(6, plotOutput("scatter_rain_snow", height = "340px", click = "scatter_click")),
        column(6, plotOutput("cdf_click", height = "340px"))
      )
    )
  )
)

server <- function(input, output, session) {
  peaks <- reactiveVal(NULL)
  dl_models <- reactiveVal(NULL)
  pp_tbl <- reactiveVal(NULL)
  skill_tbl <- reactiveVal(NULL)
  map_summary <- reactiveVal(NULL)
  diag_error <- reactiveVal(NULL)
  diag_source <- reactiveVal(NA_character_)
  cdf_click <- reactiveVal(NULL)

  output$data_banner <- renderUI({
    missing <- c()
    if (!peaks_ok) {
      missing <- c(missing, basename(peaks_path))
    }
    if (!models_ok) {
      missing <- c(missing, basename(models_path))
    }
    if (length(missing) == 0L) {
      return(NULL)
    }
    div(
      class = "alert alert-warning",
      strong("Missing: "),
      paste(missing, collapse = ", "),
      " — run script 4 from ",
      code("apps/4a-distributional-learning-fit"),
      " first."
    )
  })

  load_models_into_session <- function() {
    if (!models_ok) {
      return(FALSE)
    }
    tryCatch(
      {
        dl_models(readRDS(models_path))
        TRUE
      },
      error = function(e) {
        showNotification(paste("Could not read models:", conditionMessage(e)), type = "error", duration = NULL)
        FALSE
      }
    )
  }

  observe({
    if (!peaks_ok) {
      return(invisible(NULL))
    }
    if (is.null(peaks())) {
      peaks(readRDS(peaks_path))
      return(invisible(NULL))
    }
    if (is.null(dl_models()) && models_ok) {
      load_models_into_session()
    }
  })

  output$load_status <- renderText({
    lines <- character()
    if (!peaks_ok) {
      return("Waiting for POT peaks (script 3).")
    }
    if (!models_ok) {
      return("Waiting for fitted models (run script 4).")
    }
    if (!is.null(dl_models())) {
      lines <- c(lines, paste0("Models loaded (", nrow(dl_models()), " cells). CDF clicks ready."))
    } else {
      lines <- c(lines, "Models not loaded yet.")
    }
    if (!is.null(pp_tbl())) {
      src <- switch(
        diag_source() %||% "",
        saved = "from diagnostics RDS",
        recomputed = "recomputed from predictions RDS",
        rebuilt = "rebuilt from models + peaks",
        "loaded"
      )
      lines <- c(lines, paste0("P–P / skill / map loaded (", src, ")."))
    } else {
      err <- diag_error()
      if (length(err) && nzchar(err)) {
        lines <- c(lines, paste("Diagnostics error:", err))
      } else {
        lines <- c(lines, "Click Load P–P and skill diagnostics for calibration plots and map.")
      }
    }
    paste(lines, collapse = "\n")
  })

  observeEvent(input$reload_models, {
    if (load_models_into_session()) {
      showNotification("Models reloaded from disk.", type = "message")
    }
  })

  observeEvent(input$load_diagnostics, {
    req(peaks())
    diag_error(NULL)

    tryCatch(
      {
        withProgress(message = "Loading diagnostics…", value = 0, {
          incProgress(0.2, detail = if (diagnostics_ok) "Diagnostics RDS" else "Computing")
          res <- dl_load_diagnostics(
            repo_root,
            peaks(),
            dl_models = dl_models(),
            diagnostics_path = diag_path,
            predictions_path = dl_path
          )
          if (!isTRUE(res$ok)) {
            stop(res$message %||% "Could not load diagnostics.", call. = FALSE)
          }
          if (identical(res$source, "rebuilt") && nzchar(res$predictions_error %||% "")) {
            showNotification(
              paste(
                "Could not read", basename(dl_path), "(", res$predictions_error, ").",
                "Recomputed diagnostics from models and peaks."
              ),
              type = "warning",
              duration = 10
            )
          } else if (identical(res$source, "recomputed")) {
            showNotification(
              paste(
                "No diagnostics RDS — recomputed from",
                basename(dl_path), "."
              ),
              type = "warning",
              duration = 8
            )
          }
          incProgress(1, detail = "Done")
          pp_tbl(res$pp)
          skill_tbl(res$skill)
          map_summary(res$map_summary)
          diag_source(res$source)
        })
      },
      error = function(e) {
        diag_error(conditionMessage(e))
        showNotification(paste("Diagnostics failed:", conditionMessage(e)), type = "error", duration = NULL)
      }
    )
  })

  cells_ref <- reactive({
    req(peaks())
    peaks() |>
      dplyr::distinct(cell_id, x, y) |>
      dplyr::arrange(cell_id)
  })

  observe({
    cr <- cells_ref()
    updateSelectInput(session, "cell_id", choices = cr$cell_id, selected = min(cr$cell_id))
  })

  observeEvent(input$cell_id, {
    cdf_click(NULL)
  })

  selected_cell_id <- reactive({
    req(input$cell_id)
    as.integer(input$cell_id)
  })

  output$cell_meta <- renderText({
    req(peaks(), input$cell_id)
    row <- cells_ref() |> dplyr::filter(.data$cell_id == as.integer(input$cell_id))
    paste0(
      "lon = ", round(row$y, 3), ", lat = ", round(row$x, 3),
      "\nModels: ", if (!is.null(dl_models())) "loaded" else "missing",
      " | Diagnostics: ", if (!is.null(pp_tbl())) "loaded" else "not loaded"
    )
  })

  map_xlim <- reactive(range(cells_ref()$y, na.rm = TRUE) + c(-0.5, 0.5))
  map_ylim <- reactive(range(cells_ref()$x, na.rm = TRUE) + c(-0.5, 0.5))

  world_map <- reactive({
    bb <- st_bbox(
      c(xmin = map_xlim()[1], xmax = map_xlim()[2], ymin = map_ylim()[1], ymax = map_ylim()[2]),
      crs = st_crs(4326)
    )
    rnaturalearth::ne_countries(scale = 50, returnclass = "sf") |> st_crop(bb)
  })

  td <- reactive(dl_tile_dims(cells_ref()))

  observeEvent(input$map_click, {
    req(cells_ref())
    cid <- dl_nearest_cell(input$map_click$x, input$map_click$y, cells_ref())
    updateSelectInput(session, "cell_id", selected = cid)
  })

  output$map_skill <- renderPlot({
    req(!is.null(map_summary()))
    sel <- selected_cell_id()
    d <- map_summary() |> dplyr::mutate(selected = .data$cell_id == sel)
    ggplot(d, aes(y, x)) +
      geom_sf(data = world_map(), inherit.aes = FALSE, fill = NA, linewidth = 1) +
      geom_tile(
        aes(fill = skill_median),
        linewidth = 0.25,
        colour = scales::alpha("grey35", 0.35),
        width = td()["width"],
        height = td()["height"],
        alpha = 0.88
      ) +
      geom_tile(
        data = dplyr::filter(d, selected),
        aes(y, x),
        inherit.aes = FALSE,
        fill = NA,
        colour = "grey10",
        linewidth = 0.65,
        width = td()["width"],
        height = td()["height"]
      ) +
      geom_text(aes(label = cell_id), size = 2.2, colour = "grey15") +
      scale_fill_gradientn(
        "Median skill\n(vs marginal)",
        colours = dl_pal,
        na.value = "grey92",
        labels = scales::percent_format()
      ) +
      coord_sf(xlim = map_xlim(), ylim = map_ylim(), expand = FALSE) +
      labs(
        x = "Longitude", y = "Latitude",
        title = "Median quantile skill (all τ)"
      ) +
      theme_minimal() +
      theme(panel.grid = element_blank())
  })

  output$pp_plot <- renderPlot({
    req(!is.null(pp_tbl()))
    sel <- selected_cell_id()
    d <- pp_tbl()
    ggplot(aes(p_empirical, p_model)) +
      geom_line(
        data = dplyr::filter(d, .data$cell_id != sel),
        aes(group = interaction(cell_id, model)),
        alpha = 0.12,
        linewidth = 0.35
      ) +
      geom_line(
        data = dplyr::filter(d, .data$cell_id == sel),
        linewidth = 0.9,
        colour = "firebrick"
      ) +
      geom_abline(intercept = 0, slope = 1, linetype = "dashed", colour = "orange3") +
      facet_wrap(~model, nrow = 1, labeller = labeller(
        model = c(forest = "Forest (empirical)", gp = "GP tail")
      )) +
      labs(title = paste0("P–P calibration | cell ", sel)) +
      theme_bw()
  })

  output$skill_plot <- renderPlot({
    req(!is.null(skill_tbl()))
    sel <- selected_cell_id()
    d <- skill_tbl()
    ggplot(aes(tau, skill_score)) +
      geom_line(
        data = dplyr::filter(d, .data$cell_id != sel),
        aes(group = cell_id),
        alpha = 0.18,
        linewidth = 0.35
      ) +
      geom_line(
        data = dplyr::filter(d, .data$cell_id == sel),
        colour = "firebrick",
        linewidth = 0.9
      ) +
      scale_y_continuous("Skill", labels = scales::percent_format()) +
      labs(title = paste0("Skill vs marginal | cell ", sel), x = "τ") +
      theme_bw()
  })

  peaks_cell <- reactive({
    req(peaks(), input$cell_id)
    peaks() |>
      dplyr::filter(.data$cell_id == as.integer(input$cell_id)) |>
      dplyr::select(rainfall_hourly, snowmelt_hourly, runoff_hourly, date)
  })

  output$scatter_rain_snow <- renderPlot({
    req(nrow(peaks_cell()) > 0)
    ggplot(peaks_cell(), aes(rainfall_hourly, snowmelt_hourly)) +
      geom_point(alpha = 0.35, size = 1.1, colour = "grey25") +
      coord_cartesian(expand = FALSE) +
      labs(
        title = paste0("Rain vs snowmelt (peak hours) | cell ", input$cell_id),
        x = "Rainfall (mm/h)", y = "Snowmelt (mm/h)"
      ) +
      theme_bw()
  })

  observeEvent(input$scatter_click, {
    req(input$cell_id, nrow(peaks_cell()) > 0)
    model <- dl_get_cell_model(as.integer(input$cell_id), dl_models())
    if (is.null(model)) {
      showNotification(
        "No model for this cell — run script 4 or reload models.",
        type = "warning",
        duration = 6
      )
      return(invisible(NULL))
    }

    clk <- input$scatter_click
    df <- peaks_cell()
    dxy <- sqrt((df$rainfall_hourly - clk$x)^2 + (df$snowmelt_hourly - clk$y)^2)
    min_i <- which.min(dxy)
    span <- sqrt(diff(range(df$rainfall_hourly))^2 + diff(range(df$snowmelt_hourly))^2)
    runoff_obs <- if (dxy[[min_i]] <= max(1e-9, 0.02 * span)) df$runoff_hourly[[min_i]] else NA_real_

    dst_f <- predict(model, newdata = tibble(rainfall_hourly = clk$x, snowmelt_hourly = clk$y))[[1]]
    dst_g <- tryCatch(convert_emp_to_gp(dst_f), error = function(e) dst_f)
    br <- range(df$runoff_hourly, na.rm = TRUE)
    grid <- seq(max(0, br[1] * 0.85), max(br[2] * 1.12, br[1] + 1e-6), length.out = 400)

    cdf_click(list(
      rain = clk$x,
      snow = clk$y,
      runoff_obs = runoff_obs,
      grid = grid,
      cdf_forest = as.numeric(distionary::eval_cdf(dst_f, at = grid)),
      cdf_gp = as.numeric(distionary::eval_cdf(dst_g, at = grid))
    ))
  })

  output$cdf_click <- renderPlot({
    z <- cdf_click()
    if (is.null(z)) {
      return(ggplot() +
        annotate("text", x = 0.5, y = 0.55, label = "Click the scatter plot at a rainfall–snowmelt pair.") +
        theme_void())
    }
    cdf_tbl <- tibble(runoff = z$grid, forest = z$cdf_forest, gp = z$cdf_gp) |>
      tidyr::pivot_longer(-runoff, names_to = "version", values_to = "cdf") |>
      dplyr::mutate(version = factor(version, c("forest", "gp"),
        labels = c("Random forest (empirical)", "GP tail")))
    p <- ggplot(cdf_tbl, aes(runoff, cdf, colour = version)) +
      geom_line(linewidth = 0.85) +
      scale_colour_manual(values = c("grey25", "steelblue")) +
      coord_cartesian(xlim = range(z$grid), ylim = c(0, 1), expand = FALSE) +
      labs(
        title = sprintf("Conditional runoff CDF | rain = %.3g, snow = %.3g mm/h", z$rain, z$snow),
        x = "Hourly runoff (mm)", y = "CDF", colour = NULL
      ) +
      theme_bw() + theme(legend.position = "bottom")
    if (is.finite(z$runoff_obs)) {
      p <- p + geom_vline(xintercept = z$runoff_obs, linetype = "dashed", colour = "red3")
    }
    p
  })
}

shinyApp(ui, server)

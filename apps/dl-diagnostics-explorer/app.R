# Distributional-learning diagnostics: calibration (P–P), skill vs marginal,
# rain–snow scatter with conditional runoff CDFs at clicked coordinates.
# Run from the package root:
#   shiny::runApp("apps/dl-diagnostics-explorer")
#
# Requires:
#   - data/era5_land_hourly_alps_peaks.rds
#   - data/era5_land_hourly_alps_dl_predictions.rds  (script 4)
#   - data/era5_land_hourly_alps_dl_rqforest_models.rds  (script 4)

library(shiny)
library(tidyverse)
library(sf)
library(fs)
library(rnaturalearth)
library(probaverse)
library(distionary)

repo_root <- here::here()
devtools::load_all(repo_root, quiet = TRUE)

peaks_path <- path(repo_root, "data", "era5_land_hourly_alps_peaks.rds")
dl_path <- path(repo_root, "data", "era5_land_hourly_alps_dl_predictions.rds")
models_path <- path(repo_root, "data", "era5_land_hourly_alps_dl_rqforest_models.rds")

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

compute_pp_long <- function(peak_hour_distributions) {
  peak_hour_distributions |>
    group_by(cell_id, x, y) |>
    mutate(
      p_model_forest = map2_dbl(distribution_forest, runoff_hourly, eval_cdf),
      p_empirical_forest = uscore(p_model_forest),
      p_model_gp = map2_dbl(distribution_gp, runoff_hourly, eval_cdf),
      p_empirical_gp = uscore(p_model_gp)
    ) |>
    ungroup() |>
    select(cell_id, x, y, starts_with("p_")) |>
    pivot_longer(
      starts_with("p_"),
      names_to = c(".value", "model"),
      names_pattern = "(p_.*)_(.*)"
    ) |>
    arrange(cell_id, model, p_empirical)
}

compute_skill_scores <- function(peak_hour_distributions, dat) {
  qscores_model <- peak_hour_distributions |>
    mutate(
      df = map(
        distribution_forest,
        enframe_quantile,
        at = 1:99 / 100,
        arg_name = "tau"
      )
    ) |>
    unnest(df) |>
    group_by(cell_id, x, y, tau) |>
    summarise(
      qscore_model = mean(quantile_score(
        runoff_hourly,
        xhat = quantile,
        tau = tau
      )),
      .groups = "drop"
    )

  null_model <- dat |>
    group_by(cell_id, x, y) |>
    summarise(runoff_hourly = list(runoff_hourly), .groups = "drop") |>
    mutate(marginal = map(runoff_hourly, dst_empirical))

  null_quantiles <- null_model |>
    mutate(
      df = map(marginal, enframe_quantile, at = 1:99 / 100, arg_name = "tau")
    ) |>
    select(!marginal) |>
    unnest(df)

  qscores_null <- null_quantiles |>
    mutate(
      qscore_null = pmap_dbl(
        list(runoff_hourly, quantile, tau),
        \(y, q, p) mean(quantile_score(y, xhat = q, tau = p))
      )
    ) |>
    select(x, y, tau, qscore_null)

  left_join(qscores_null, qscores_model, by = c("x", "y", "tau")) |>
    mutate(skill_score = 1 - qscore_model / qscore_null)
}

data_ready <- file.exists(peaks_path) &&
  file.exists(dl_path) &&
  file.exists(models_path)

ui <- fluidPage(
  titlePanel("Distributional learning diagnostics"),
  uiOutput("data_banner"),
  helpText(
    "Inspect P–P calibration and quantile skill (vs the cell’s marginal empirical forecast) ",
    "for each grid cell. Click the map to pick a cell. In the rain–snow scatter plot, ",
    "click anywhere to view the predicted hourly runoff CDF at that rainfall and snowmelt rate; ",
    "if your click lands on an observed peak (or very close), the realized runoff is shown ",
    "as a vertical dashed line. There are no tuning inputs yet (future: hyperparameters)."
  ),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      selectInput(
        "cell_id",
        "Cell ID",
        choices = integer(0),
        selected = NULL
      ),
      verbatimTextOutput("cell_meta"),
      verbatimTextOutput("prep_status")
    ),
    mainPanel(
      width = 9,
      plotOutput("map_skill", height = "320px", click = "map_click"),
      fluidRow(
        column(6, plotOutput("pp_plot", height = "300px")),
        column(6, plotOutput("skill_plot", height = "300px"))
      ),
      fluidRow(
        column(
          6,
          plotOutput("scatter_rain_snow", height = "340px", click = "scatter_click")
        ),
        column(6, plotOutput("cdf_click", height = "340px"))
      )
    )
  )
)

server <- function(input, output, session) {
  values <- reactiveValues(
    ready = FALSE,
    pp_tbl = NULL,
    skill_tbl = NULL,
    map_summary = NULL,
    peaks = NULL,
    dl_models = NULL,
    dl_predictions = NULL,
    error = NULL,
    cdf_click = NULL
  )

  output$data_banner <- renderUI({
    if (data_ready) {
      return(NULL)
    }
    missing <- c()
    if (!file.exists(peaks_path)) {
      missing <- c(missing, basename(peaks_path))
    }
    if (!file.exists(dl_path)) {
      missing <- c(missing, basename(dl_path))
    }
    if (!file.exists(models_path)) {
      missing <- c(missing, basename(models_path))
    }
    helpText(
      strong("Missing data: "),
      paste(missing, collapse = "; "),
      ". Run ",
      code("scripts/4-distributional_learning.r"),
      ", then reload."
    )
  })

  observe({
    if (!data_ready) {
      return(invisible(NULL))
    }

    values$error <- NULL
    values$ready <- FALSE
    values$cdf_click <- NULL

    tryCatch(
      {
        withProgress(
          message = "Loading distributional learning outputs…",
          detail = "Reading RDS",
          value = 0,
          {
            incProgress(0.1, detail = "Peaks")
            peaks <- read_rds(peaks_path)

            incProgress(0.25, detail = "Hourly DL predictions")
            dl_predictions <- read_rds(dl_path)

            incProgress(0.35, detail = "Fitted rqforest models")
            dl_models <- read_rds(models_path)

            incProgress(0.45, detail = "P–P table")
            pp_tbl <- compute_pp_long(dl_predictions)

            incProgress(0.65, detail = "Skill scores")
            skill_tbl <- compute_skill_scores(dl_predictions, peaks)

            incProgress(0.82, detail = "Map summary")
            map_summary <- skill_tbl |>
              group_by(cell_id, x, y) |>
              summarise(
                skill_median = median(skill_score, na.rm = TRUE),
                .groups = "drop"
              )

            incProgress(1, detail = "Done")
            values$peaks <- peaks
            values$dl_predictions <- dl_predictions
            values$dl_models <- dl_models
            values$pp_tbl <- pp_tbl
            values$skill_tbl <- skill_tbl
            values$map_summary <- map_summary
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

  cells_ref <- reactive({
    req(values$ready, values$peaks)
    values$peaks |>
      distinct(cell_id, x, y) |>
      arrange(cell_id)
  })

  observe({
    cr <- cells_ref()
    updateSelectInput(
      session,
      "cell_id",
      choices = cr$cell_id,
      selected = min(cr$cell_id)
    )
  })

  observeEvent(input$cell_id, {
    values$cdf_click <- NULL
  })

  map_xlim <- reactive({
    cr <- cells_ref()
    range(cr$y, na.rm = TRUE) + c(-0.5, 0.5)
  })

  map_ylim <- reactive({
    cr <- cells_ref()
    range(cr$x, na.rm = TRUE) + c(-0.5, 0.5)
  })

  map_bbox <- reactive({
    req(map_xlim(), map_ylim())
    st_bbox(
      c(
        xmin = map_xlim()[1],
        xmax = map_xlim()[2],
        ymin = map_ylim()[1],
        ymax = map_ylim()[2]
      ),
      crs = st_crs(4326)
    )
  })

  world_map <- reactive({
    req(map_bbox())
    ne_countries(scale = 50, returnclass = "sf") |>
      st_crop(map_bbox())
  })

  td <- reactive({
    req(cells_ref())
    tile_dims(cells_ref())
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

  output$cell_meta <- renderText({
    req(values$ready, input$cell_id)
    row <- cells_ref() |> filter(cell_id == as.integer(input$cell_id))
    if (nrow(row) != 1L) {
      return("")
    }
    paste0(
      "lon = ", round(row$y, 3), ", lat = ", round(row$x, 3),
      "\nPOT peak hours per cell; scatter shows rainfall vs snowmelt for peaks."
    )
  })

  observeEvent(input$map_click, {
    req(cells_ref())
    cid <- nearest_cell(input$map_click$x, input$map_click$y, cells_ref())
    updateSelectInput(session, "cell_id", selected = cid)
  })

  output$map_skill <- renderPlot({
    req(values$ready, values$map_summary)
    d <- values$map_summary |>
      mutate(selected = cell_id == as.integer(input$cell_id))

    ggplot(d, aes(y, x)) +
      geom_sf(
        data = world_map(),
        inherit.aes = FALSE,
        fill = NA,
        linewidth = 1
      ) +
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
        colours = pal,
        na.value = "grey92",
        labels = scales::percent_format()
      ) +
      coord_sf(xlim = map_xlim(), ylim = map_ylim(), expand = FALSE) +
      labs(
        x = "Longitude",
        y = "Latitude",
        title = "Median quantile skill score (all τ), marginal empirical reference"
      ) +
      theme_minimal() +
      theme(panel.grid = element_blank())
  })

  output$pp_plot <- renderPlot({
    req(values$ready, values$pp_tbl, input$cell_id)
    sel <- as.integer(input$cell_id)
    d <- values$pp_tbl
    d_bg <- d |> filter(cell_id != sel)
    d_hi <- d |> filter(cell_id == sel)

    ggplot(mapping = aes(p_empirical, p_model)) +
      geom_line(
        data = d_bg,
        aes(group = interaction(cell_id, model)),
        alpha = 0.12,
        linewidth = 0.35
      ) +
      geom_line(
        data = d_hi,
        linewidth = 0.9,
        alpha = 0.95,
        colour = "firebrick"
      ) +
      geom_abline(
        intercept = 0,
        slope = 1,
        linetype = "dashed",
        colour = "orange3"
      ) +
      facet_wrap(~model, nrow = 1, labeller = labeller(
        model = c(forest = "Forest (empirical)", gp = "GP tail")
      )) +
      labs(
        x = expression(hat(u) ~ "(empirical among hours)"),
        y = expression(F(y) ~ "from forecast"),
        title = paste0("P–P calibration | cell ", sel)
      ) +
      theme_bw() +
      theme(legend.position = "none")
  })

  output$skill_plot <- renderPlot({
    req(values$ready, values$skill_tbl, input$cell_id)
    sel <- as.integer(input$cell_id)
    d <- values$skill_tbl
    d_bg <- d |> filter(cell_id != sel)
    d_hi <- d |> filter(cell_id == sel)

    ggplot(mapping = aes(tau, skill_score)) +
      geom_line(
        data = d_bg,
        aes(group = cell_id),
        alpha = 0.18,
        linewidth = 0.35
      ) +
      geom_line(
        data = d_hi,
        colour = "firebrick",
        linewidth = 0.9
      ) +
      labs(
        x = "Quantile level τ",
        title = paste0("Skill vs marginal | cell ", sel)
      ) +
      scale_y_continuous("Skill", labels = scales::percent_format()) +
      theme_bw()
  })

  peaks_cell <- reactive({
    req(values$ready, values$peaks, input$cell_id)
    values$peaks |>
      filter(cell_id == as.integer(input$cell_id)) |>
      select(rainfall_hourly, snowmelt_hourly, runoff_hourly, date)
  })

  output$scatter_rain_snow <- renderPlot({
    req(nrow(peaks_cell()) > 0)
    df <- peaks_cell()
    ggplot(df, aes(rainfall_hourly, snowmelt_hourly)) +
      geom_point(alpha = 0.35, size = 1.1, colour = "grey25") +
      coord_cartesian(expand = FALSE) +
      labs(
        title = paste0("Rain vs snowmelt (peak hours) | cell ", input$cell_id),
        x = "Rainfall (mm/h)",
        y = "Snowmelt (mm/h)"
      ) +
      theme_bw()
  })

  observeEvent(input$scatter_click, {
    req(values$ready, input$cell_id, nrow(peaks_cell()) > 0)

    clk <- input$scatter_click
    rain_click <- clk$x
    snow_click <- clk$y

    df <- peaks_cell()
    rx <- diff(range(df$rainfall_hourly, na.rm = TRUE))
    ry <- diff(range(df$snowmelt_hourly, na.rm = TRUE))
    span <- sqrt(rx^2 + ry^2)
    thr <- max(1e-9, 0.02 * span)

    dxy <- sqrt((df$rainfall_hourly - rain_click)^2 + (df$snowmelt_hourly - snow_click)^2)
    min_i <- which.min(dxy)
    runoff_at_click <- if (dxy[[min_i]] <= thr) {
      df$runoff_hourly[[min_i]]
    } else {
      NA_real_
    }

    mod_tbl <- values$dl_models
    row <- mod_tbl |> dplyr::filter(cell_id == as.integer(input$cell_id))
    if (nrow(row) != 1L) {
      showNotification("No rqforest model for this cell.", type = "warning", duration = 6)
      return(invisible(NULL))
    }

    dst_f <- predict(row$dl_rqforest[[1]], newdata = tibble(
      rainfall_hourly = rain_click,
      snowmelt_hourly = snow_click
    ))[[1]]

    dst_g <- tryCatch(
      convert_emp_to_gp(dst_f),
      error = function(e) dst_f
    )

    br <- range(df$runoff_hourly, na.rm = TRUE)
    grid <- seq(max(0, br[1] * 0.85), max(br[2] * 1.12, br[1] + 1e-6), length.out = 400)

    cdf_forest <- tryCatch(
      distionary::eval_cdf(dst_f, at = grid),
      error = function(e) {
        showNotification(
          paste("Forest CDF evaluation failed:", conditionMessage(e)),
          type = "warning",
          duration = 8
        )
        rep(NA_real_, length(grid))
      }
    )
    cdf_gp <- tryCatch(
      distionary::eval_cdf(dst_g, at = grid),
      error = function(e) {
        showNotification(
          paste("GP CDF evaluation failed:", conditionMessage(e)),
          type = "warning",
          duration = 8
        )
        rep(NA_real_, length(grid))
      }
    )

    values$cdf_click <- list(
      rain = rain_click,
      snow = snow_click,
      runoff_obs = runoff_at_click,
      grid = grid,
      cdf_forest = cdf_forest,
      cdf_gp = cdf_gp
    )
  })

  output$cdf_click <- renderPlot({
    if (is.null(values$cdf_click)) {
      return(
        ggplot() +
          annotate(
            "text",
            x = 0.5,
            y = 0.55,
            label = "Click the scatter plot at a rainfall–snowmelt pair."
          ) +
          coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
          theme_void()
      )
    }
    z <- values$cdf_click
    cdf_tbl <- tibble(
      runoff = z$grid,
      forest = z$cdf_forest,
      gp = z$cdf_gp
    ) |>
      pivot_longer(-runoff, names_to = "version", values_to = "cdf") |>
      mutate(
        version = factor(
          version,
          levels = c("forest", "gp"),
          labels = c("Random forest (empirical)", "GP tail")
        )
      )

    p <- ggplot(cdf_tbl, aes(runoff, cdf, colour = version)) +
      geom_line(linewidth = 0.85) +
      scale_colour_manual(values = c("grey25", "steelblue")) +
      coord_cartesian(xlim = range(z$grid), ylim = c(0, 1), expand = FALSE) +
      labs(
        x = "Hourly runoff (mm)",
        y = "CDF",
        title = sprintf(
          "Conditional runoff CDF | rain = %.3g, snow = %.3g mm/h",
          z$rain,
          z$snow
        ),
        colour = NULL
      ) +
      theme_bw() +
      theme(legend.position = "bottom")

    if (is.finite(z$runoff_obs)) {
      p <- p +
        geom_vline(
          xintercept = z$runoff_obs,
          linetype = "dashed",
          colour = "red3",
          linewidth = 0.45
        ) +
        labs(subtitle = sprintf("Dashed: realized runoff %.4g mm/h (nearest peak)", z$runoff_obs))
    } else {
      p <- p + labs(subtitle = "Click closer to a dot to mark observed runoff at that peak.")
    }
    p
  })
}

shinyApp(ui, server)

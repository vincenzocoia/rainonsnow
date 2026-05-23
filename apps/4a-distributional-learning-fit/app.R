# Distributional learning — fit / tune (script 4).
# Run from the package root:
#   shiny::runApp("apps/4a-distributional-learning-fit")
#
# Requires: derived/era5_land_hourly_alps_peaks.rds (script 3).
# Writes: inputs/distributional_learning.yaml; runs scripts/4-distributional_learning.r
# for all cells. Inspect one cell in-app, then open 4b for full-grid diagnostics.

library(shiny)
library(tidyverse)
library(yaml)
library(probaverse)
library(distionary)

repo_root <- here::here()
devtools::load_all(repo_root, quiet = TRUE)
source(fs::path(repo_root, "apps", "dl_shared.R"))

meta_path <- fs::path(repo_root, "inputs", "distributional_learning.yaml")
peaks_path <- fs::path(repo_root, "derived", "era5_land_hourly_alps_peaks.rds")

dl0 <- dl_read_cfg(meta_path)
peaks_ok <- file.exists(peaks_path)
na_action_choices <- c("drop", "null")

empty_preview_plot <- function(message) {
  ggplot2::ggplot() +
    ggplot2::annotate(
      "text",
      x = 0.5,
      y = 0.55,
      label = message,
      size = 4.2,
      lineheight = 0.95
    ) +
    ggplot2::coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
    ggplot2::theme_void()
}

ui <- fluidPage(
  tags$head(
    tags$style(HTML("
      @keyframes dl-fit-spin { to { transform: rotate(360deg); } }
      .dl-fit-spinner {
        display: inline-block;
        animation: dl-fit-spin 0.9s linear infinite;
        margin-right: 6px;
      }
    "))
  ),
  titlePanel("Distributional learning — fit (script 4a)"),
  uiOutput("data_banner"),
  helpText(
    "Tune ",
    code("dl_rqforest"),
    " and preview on ",
    strong("one cell"),
    ". When satisfied, ",
    strong("Run script 4 for all cells"),
    " writes ",
    code("inputs/distributional_learning.yaml"),
    " and derived outputs. Open ",
    code("apps/4b-distributional-learning-diagnostics"),
    " to explore the fitted grid."
  ),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      tags$p(
        style = "font-size: 0.9em; margin-bottom: 14px;",
        strong("Model (fixed): "),
        "hourly runoff ~ rainfall + snowmelt on POT peak hours."
      ),
      selectInput("na_action", "Missing values", na_action_choices, dl0$na_action %||% "drop"),
      dl_param_hint(tags$span(
        tags$strong("drop"), " removes rows with missing values; ",
        tags$strong("null"), " uses a null model if any NA exists in the cell."
      )),
      numericInput("min_obs", "Minimum observations per cell", dl0$min_obs %||% 5L, 1L, step = 1L),
      dl_param_hint(tags$span(
        tags$strong("Lower"), " fits more cells; ",
        tags$strong("higher"), " is stricter and stabler per cell."
      )),
      numericInput("ntree", "Number of trees (ntree)", dl0$ntree %||% 500L, 10L, step = 10L),
      dl_param_hint(tags$span(
        tags$strong("More trees"), " = smoother, slower; ",
        tags$strong("fewer"), " = faster, noisier."
      )),
      numericInput("mtry", "Predictors per split (mtry)", dl0$mtry %||% 1L, 1L, max = 2L, step = 1L),
      dl_param_hint(tags$span(
        tags$strong("mtry = 1"), " is simpler; ",
        tags$strong("mtry = 2"), " allows rain and snow together at splits."
      )),
      numericInput("nodesize", "Minimum leaf size (nodesize)", dl0$nodesize %||% 5L, 1L, step = 1L),
      dl_param_hint(tags$span(
        tags$strong("Smaller"), " leaves = more flexible; ",
        tags$strong("larger"), " leaves = simpler surfaces."
      )),
      numericInput("nthreads", "Parallel threads (nthreads)", dl0$nthreads %||% 1L, 1L, step = 1L),
      dl_param_hint(tags$span(tags$strong("More threads"), " speeds up fitting only.")),
      hr(),
      selectInput("cell_id", "Cell ID (preview)", choices = integer(0)),
      actionButton("fit_cell", "Fit selected cell (preview)", class = "btn-info", width = "100%"),
      uiOutput("preview_status_ui"),
      actionButton("run_script4", "Run script 4 for all cells", class = "btn-primary", width = "100%"),
      verbatimTextOutput("run_log"),
      verbatimTextOutput("cell_meta")
    ),
    mainPanel(
      width = 9,
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
  preview <- reactiveVal(NULL)
  cdf_click <- reactiveVal(NULL)
  fit_busy <- reactiveVal(FALSE)
  fit_error <- reactiveVal(NULL)
  run_log_txt <- reactiveVal("Adjust hyperparameters, preview one cell, then run script 4.")

  output$data_banner <- renderUI({
    if (peaks_ok) {
      return(NULL)
    }
    div(
      class = "alert alert-warning",
      strong("Missing "),
      basename(peaks_path),
      " — run script 3 (POT) first."
    )
  })

  observe({
    if (!peaks_ok || !is.null(peaks())) {
      return(invisible(NULL))
    }
    peaks(readRDS(peaks_path))
  })

  observe({
    cfg <- dl_read_cfg(meta_path)
    updateSelectInput(session, "na_action", selected = cfg$na_action %||% "drop")
    updateNumericInput(session, "min_obs", value = cfg$min_obs %||% 5L)
    updateNumericInput(session, "ntree", value = cfg$ntree %||% 500L)
    updateNumericInput(session, "mtry", value = cfg$mtry %||% 1L)
    updateNumericInput(session, "nodesize", value = cfg$nodesize %||% 5L)
    updateNumericInput(session, "nthreads", value = cfg$nthreads %||% 1L)
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
    fit_error(NULL)
  })

  output$cell_meta <- renderText({
    req(peaks(), input$cell_id)
    row <- cells_ref() |> dplyr::filter(.data$cell_id == as.integer(input$cell_id))
    paste0(
      "lon = ", round(row$y, 3), ", lat = ", round(row$x, 3),
      "\nPOT peak hours for preview."
    )
  })

  observeEvent(input$fit_cell, {
    if (isTRUE(fit_busy())) {
      return(invisible(NULL))
    }
    req(peaks(), input$cell_id)

    cid <- as.integer(input$cell_id)
    fit_busy(TRUE)
    fit_error(NULL)
    preview(NULL)
    cdf_click(NULL)
    updateActionButton(session, "fit_cell", label = "Fitting…")

    on.exit({
      fit_busy(FALSE)
      updateActionButton(session, "fit_cell", label = "Fit selected cell (preview)")
    }, add = TRUE)

    res <- NULL
    t0 <- Sys.time()
    tryCatch(
      {
        withProgress(
          message = paste0("Fitting preview | cell ", cid),
          value = 0,
          min = 0,
          max = 1,
          {
            incProgress(0.15, detail = "Training rqforest")
            res <- dl_fit_cell_preview(peaks(), cid, dl_cfg_from_input(input))
            incProgress(0.85, detail = "P–P and skill tables")
            incProgress(1, detail = "Done")
          }
        )
      },
      error = function(e) {
        fit_error(conditionMessage(e))
        showNotification(
          paste("Preview fit failed:", conditionMessage(e)),
          type = "error",
          duration = NULL
        )
      }
    )

    elapsed <- round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1)

    if (is.null(res) || !isTRUE(res$ok)) {
      msg <- if (!is.null(res) && nzchar(res$message %||% "")) {
        res$message
      } else {
        fit_error() %||% "Preview fit did not complete."
      }
      fit_error(msg)
      showNotification(msg, type = "warning", duration = 8)
      return(invisible(NULL))
    }

    res$fit_seconds <- elapsed
    preview(res)
    cdf_click(NULL)

    showNotification(
      tags$div(
        tags$strong("Preview fit complete"),
        tags$br(),
        sprintf(
          "Cell %d — %d peak hours, median skill %s (%.1f s). Plots updated.",
          res$cell_id,
          res$n_peaks,
          sprintf("%.1f%%", 100 * res$skill_median),
          elapsed
        )
      ),
      type = "message",
      duration = 10
    )
  })

  output$preview_status_ui <- renderUI({
    if (isTRUE(fit_busy())) {
      return(div(
        class = "alert alert-info",
        style = "padding: 10px; margin-top: 8px;",
        tags$span(class = "dl-fit-spinner", "\u27f3"),
        tags$strong(paste0("Fitting cell ", input$cell_id, "…")),
        tags$br(),
        tags$span(
          style = "font-size: 0.9em;",
          "Training the forest and building preview plots (often 10–60 s with many trees)."
        )
      ))
    }

    err <- fit_error()
    if (length(err) && nzchar(err)) {
      return(div(
        class = "alert alert-danger",
        style = "padding: 10px; margin-top: 8px;",
        tags$strong("Last fit failed: "),
        err
      ))
    }

    pr <- preview()
    sel <- as.integer(input$cell_id %||% NA_integer_)

    if (is.null(pr)) {
      return(div(
        class = "alert alert-secondary",
        style = "padding: 10px; margin-top: 8px;",
        "No preview yet. Choose a cell and click ",
        tags$strong("Fit selected cell"),
        " — the four plots on the right will update when finished."
      ))
    }

    stale <- !is.na(sel) && pr$cell_id != sel
    cls <- if (stale) "alert alert-warning" else "alert alert-success"
    div(
      class = cls,
      style = "padding: 10px; margin-top: 8px;",
      if (stale) {
        tagList(
          tags$strong("Preview is for a different cell."),
          tags$br(),
          sprintf(
            "Showing results for cell %d; sidebar is on cell %d — click Fit again to refresh.",
            pr$cell_id,
            sel
          )
        )
      } else {
        tagList(
          tags$strong("\u2713 Preview ready"),
          tags$br(),
          sprintf(
            "Cell %d — %d peak hours, median skill %s",
            pr$cell_id,
            pr$n_peaks,
            sprintf("%.1f%%", 100 * pr$skill_median)
          ),
          if (!is.null(pr$fit_seconds)) {
            tags$br()
            tags$span(
              style = "font-size: 0.9em;",
              sprintf(
                "Finished at %s (%.1f s). Not saved until you run script 4.",
                format(pr$fitted_at, "%H:%M:%S"),
                pr$fit_seconds
              )
            )
          } else {
            NULL
          }
        )
      }
    )
  })

  observeEvent(input$run_script4, {
    req(peaks_ok)
    dl_write_cfg(meta_path, dl_cfg_from_input(input))
    withProgress(message = "Running script 4…", value = 0.5, {
      res <- dl_run_script4(repo_root)
    })
    out_txt <- if (nzchar(res$output)) paste(res$output, res$message, sep = "\n") else res$message
    if (!isTRUE(res$ok)) {
      run_log_txt(paste0("ERROR\n", out_txt))
      showNotification(paste("Script 4 failed:", res$message), type = "error", duration = NULL)
      return(invisible(NULL))
    }
    run_log_txt(paste(
      "Last run: OK", format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
      paste0("(wrote ", basename(meta_path), " and derived/ outputs)"),
      "", out_txt, sep = "\n"
    ))
    showNotification("Script 4 finished. Open the 4b diagnostics app to explore results.", type = "message", duration = 8)
  })

  output$run_log <- renderText(run_log_txt())

  peaks_cell <- reactive({
    req(peaks(), input$cell_id)
    peaks() |>
      dplyr::filter(.data$cell_id == as.integer(input$cell_id)) |>
      dplyr::select(rainfall_hourly, snowmelt_hourly, runoff_hourly, date)
  })

  output$pp_plot <- renderPlot({
    if (isTRUE(fit_busy())) {
      return(empty_preview_plot("Fitting in progress…\nP–P plot will appear here."))
    }
    pr <- preview()
    sel <- as.integer(input$cell_id)
    if (is.null(pr) || pr$cell_id != sel) {
      return(empty_preview_plot("Click Fit selected cell\nto see P–P calibration."))
    }
    ggplot(pr$pp_tbl, aes(p_empirical, p_model)) +
      geom_line(linewidth = 0.9, colour = "firebrick") +
      geom_abline(intercept = 0, slope = 1, linetype = "dashed", colour = "orange3") +
      facet_wrap(~model, nrow = 1, labeller = labeller(
        model = c(forest = "Forest (empirical)", gp = "GP tail")
      )) +
      labs(title = paste0("P–P calibration | cell ", pr$cell_id, " (preview)")) +
      theme_bw()
  })

  output$skill_plot <- renderPlot({
    if (isTRUE(fit_busy())) {
      return(empty_preview_plot("Fitting in progress…\nSkill plot will appear here."))
    }
    pr <- preview()
    sel <- as.integer(input$cell_id)
    if (is.null(pr) || pr$cell_id != sel) {
      return(empty_preview_plot("Click Fit selected cell\nto see skill vs marginal."))
    }
    ggplot(pr$skill_tbl, aes(tau, skill_score)) +
      geom_line(colour = "firebrick", linewidth = 0.9) +
      scale_y_continuous("Skill", labels = scales::percent_format()) +
      labs(title = paste0("Skill vs marginal | cell ", pr$cell_id, " (preview)"), x = "τ") +
      theme_bw()
  })

  output$scatter_rain_snow <- renderPlot({
    req(nrow(peaks_cell()) > 0)
    ggplot(peaks_cell(), aes(rainfall_hourly, snowmelt_hourly)) +
      geom_point(alpha = 0.35, size = 1.1, colour = "grey25") +
      coord_cartesian(expand = FALSE) +
      labs(
        title = paste0("Rain vs snowmelt | cell ", input$cell_id),
        x = "Rainfall (mm/h)", y = "Snowmelt (mm/h)"
      ) +
      theme_bw()
  })

  observeEvent(input$scatter_click, {
    pr <- preview()
    req(!is.null(pr), pr$cell_id == as.integer(input$cell_id), nrow(peaks_cell()) > 0)
    clk <- input$scatter_click
    df <- peaks_cell()
    dxy <- sqrt((df$rainfall_hourly - clk$x)^2 + (df$snowmelt_hourly - clk$y)^2)
    min_i <- which.min(dxy)
    span <- sqrt(diff(range(df$rainfall_hourly))^2 + diff(range(df$snowmelt_hourly))^2)
    runoff_obs <- if (dxy[[min_i]] <= max(1e-9, 0.02 * span)) df$runoff_hourly[[min_i]] else NA_real_
    dst_f <- predict(pr$model, newdata = tibble(rainfall_hourly = clk$x, snowmelt_hourly = clk$y))[[1]]
    dst_g <- tryCatch(convert_emp_to_gp(dst_f), error = function(e) dst_f)
    br <- range(df$runoff_hourly, na.rm = TRUE)
    grid <- seq(max(0, br[1] * 0.85), max(br[2] * 1.12, br[1] + 1e-6), length.out = 400)
    cdf_click(list(
      rain = clk$x, snow = clk$y, runoff_obs = runoff_obs, grid = grid,
      cdf_forest = as.numeric(distionary::eval_cdf(dst_f, at = grid)),
      cdf_gp = as.numeric(distionary::eval_cdf(dst_g, at = grid))
    ))
  })

  output$cdf_click <- renderPlot({
    if (isTRUE(fit_busy())) {
      return(empty_preview_plot("Fitting in progress…"))
    }
    z <- cdf_click()
    pr <- preview()
    sel <- as.integer(input$cell_id)
    if (is.null(pr) || pr$cell_id != sel) {
      return(empty_preview_plot("Fit the cell, then click\nthe rain–snow scatter."))
    }
    if (is.null(z)) {
      return(empty_preview_plot("Preview ready.\nClick the scatter for a conditional CDF."))
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
        title = sprintf("Conditional runoff CDF | rain = %.3g, snow = %.3g", z$rain, z$snow),
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

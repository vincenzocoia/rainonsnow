# Rain-on-snow explorer: one app for every analysis in analyses/.
#
#   shiny::runApp("apps/explorer")
#
# Pick an analysis at the top; each tab reads that analysis's out/ folder.
# Nothing is refitted here except single predictive distributions on demand.
suppressPackageStartupMessages({
  library(shiny)
  library(bslib)
  library(tidyverse)
  library(distionary)
})
root <- normalizePath(file.path(dirname(sys.frame(1)$ofile %||% "apps/explorer/app.R"), "..", ".."),
  mustWork = FALSE)
if (!file.exists(file.path(root, "DESCRIPTION"))) root <- here::here()
devtools::load_all(root, quiet = TRUE)
source(file.path(root, "scripts", "lib", "plots.R"))

OUT_FILES <- c("training", "event_levels", "models", "diagnostics", "return_levels",
  "event_probability", "triggers", "likeliest")

load_analysis <- function(name) {
  cfg <- read_analysis(name, root = root)
  out <- list(cfg = cfg)
  for (f in OUT_FILES) {
    p <- analysis_path(name, "out", paste0(f, ".rds"), root = root)
    out[[f]] <- if (file.exists(p)) readRDS(p) else NULL
  }
  readme <- analysis_path(name, "README.md", root = root)
  out$readme <- if (file.exists(readme)) paste(readLines(readme, warn = FALSE), collapse = "\n") else NULL
  out
}

analysis_index <- function() {
  purrr::map_dfr(list_analyses(root), function(nm) {
    cfg <- read_analysis(nm, root = root)
    done <- file.exists(analysis_path(nm, "out", paste0(OUT_FILES, ".rds"), root = root))
    tibble(
      Analysis = nm,
      Title = cfg$title %||% "",
      Predictors = paste(cfg$predictors, collapse = " + "),
      Model = cfg$model$type,
      Status = cfg$status %||% "",
      Outputs = sprintf("%d / %d", sum(done), length(done))
    )
  })
}

theme <- bs_theme(version = 5, bg = "#fcfcfb", fg = "#1b2a36", primary = "#1c5cab",
  base_font = font_google("Source Sans 3"), heading_font = font_google("Source Serif 4"))

ui <- page_navbar(
  title = "Rain-on-snow explorer",
  theme = theme,
  fillable = FALSE,
  sidebar = sidebar(
    width = 290,
    selectInput("analysis", "Analysis", choices = list_analyses(root)),
    uiOutput("analysis_blurb"),
    selectInput("cell", "Grid cell", choices = NULL),
    selectInput("rp", "T-year event", choices = NULL),
    selectInput("prob", "Chance of the event", choices = c("10%" = 0.1, "50%" = 0.5, "90%" = 0.9), selected = 0.5)
  ),
  nav_panel("Overview",
    card(card_header("All analyses"), tableOutput("index")),
    layout_columns(
      card(card_header("This analysis"), uiOutput("config")),
      card(card_header("Notes (README.md)"), uiOutput("readme"))
    )
  ),
  nav_panel("Model",
    uiOutput("model_controls"),
    layout_columns(
      card(full_screen = TRUE, plotOutput("model_left", height = "480px", click = "model_click")),
      card(full_screen = TRUE, plotOutput("model_right", height = "480px"))
    )
  ),
  nav_panel("Trigger rain",
    uiOutput("trigger_controls"),
    layout_columns(
      card(full_screen = TRUE, plotOutput("trigger_left", height = "460px")),
      card(full_screen = TRUE, plotOutput("trigger_right", height = "460px"))
    ),
    card(card_header("Rain needed to trigger the event"), tableOutput("trigger_table"))
  ),
  nav_panel("Likeliest drivers",
    layout_columns(
      card(full_screen = TRUE, plotOutput("likeliest_left", height = "460px")),
      card(full_screen = TRUE, plotOutput("likeliest_right", height = "460px"))
    )
  ),
  nav_panel("Compare analyses",
    card(full_screen = TRUE, plotOutput("compare_plot", height = "480px")),
    card(card_header("Rain needed, by analysis"), tableOutput("compare_table"))
  ),
  nav_panel("Diagnostics",
    layout_columns(
      card(full_screen = TRUE, plotOutput("diag_rl", height = "420px")),
      card(full_screen = TRUE, plotOutput("diag_pp", height = "420px"))
    ),
    card(full_screen = TRUE, plotOutput("diag_skill", height = "340px"))
  )
)

server <- function(input, output, session) {
  A <- reactive({
    req(input$analysis)
    load_analysis(input$analysis)
  })
  cfg <- reactive(A()$cfg)
  p1 <- reactive(length(cfg()$predictors) == 1)
  target <- reactive(cfg()$queries$trigger$target)
  other <- reactive(setdiff(cfg()$predictors, target()))

  observeEvent(A(), {
    tr <- A()$training
    req(tr)
    cells <- sort(unique(tr$cell_id))
    sel <- if (isolate(input$cell) %in% cells) isolate(input$cell) else choose_focus_cell(cfg(), tr)
    updateSelectInput(session, "cell", choices = cells, selected = sel)
    rps <- cfg()$queries$return_periods
    sel_rp <- if (isolate(input$rp) %in% rps) isolate(input$rp) else 10
    updateSelectInput(session, "rp", choices = setNames(rps, paste0(rps, "-year")), selected = sel_rp)
  })

  cell <- reactive(as.integer(req(input$cell)))
  rp <- reactive(as.numeric(req(input$rp)))
  fc <- function(d) dplyr::filter(d, cell_id == cell())
  model <- reactive({
    m <- A()$models
    req(m)
    m$model[[match(cell(), m$cell_id)]]
  })
  cap <- reactive({
    tr <- fc(A()$training)[1, ]
    cell_label(tr$cell_id, tr$x, tr$y)
  })

  output$analysis_blurb <- renderUI({
    c <- cfg()
    tags$div(class = "small text-muted mb-2",
      tags$p(c$question),
      tags$p(tags$b("Predictors: "), paste(sapply(c$predictors, var_label), collapse = ", "),
        tags$br(), tags$b("Model: "), c$model$type))
  })

  # ---- Overview ----
  output$index <- renderTable(analysis_index())
  output$config <- renderUI({
    tags$pre(style = "font-size: 0.8rem; white-space: pre-wrap;",
      paste(readLines(analysis_path(input$analysis, "analysis.yaml", root = root)), collapse = "\n"))
  })
  output$readme <- renderUI({
    if (is.null(A()$readme)) return(tags$em("No README.md yet."))
    markdown(A()$readme)
  })

  # ---- Model ----
  output$model_controls <- renderUI({
    tr <- fc(A()$training)
    if (p1()) {
      x <- tr[[target()]]
      sliderInput("x0", paste("Show the local fit at", var_label(target())),
        min = 0, max = round(max(x), 1), value = round(unname(quantile(x, 0.75)), 1), step = 0.1, width = "100%")
    } else {
      tags$p(class = "text-muted", "Click a point on the left to see the predictive distribution at those drivers.")
    }
  })
  clicked <- reactiveVal(NULL)
  observeEvent(list(input$analysis, input$cell), clicked(NULL))
  observeEvent(input$model_click, {
    req(!p1())
    clicked(tibble(!!target() := input$model_click$x, !!other() := input$model_click$y))
  })
  output$model_left <- renderPlot({
    tr <- fc(A()$training)
    if (p1() && inherits(model(), "dl_llqr")) {
      plot_llqr_curves(model(), tr) + labs(caption = cap())
    } else {
      p <- plot_training_scatter(tr, cfg()) + facet_null() + labs(caption = cap())
      if (!is.null(clicked())) {
        p <- p + geom_point(data = clicked(), colour = ORANGE, size = 5, shape = 4, stroke = 2)
      }
      p
    }
  })
  output$model_right <- renderPlot({
    lv <- fc(A()$event_levels)
    if (p1()) {
      req(input$x0)
      if (inherits(model(), "dl_llqr")) {
        plot_llqr_local(model(), input$x0)
      } else {
        nd <- tibble(!!target() := input$x0)
        plot_predictive_exceedance(model(), nd, cfg()$tail, lv)
      }
    } else {
      tr <- fc(A()$training)
      nd <- clicked() %||% tibble(
        !!target() := unname(quantile(tr[[target()]], c(0.5, 0.9, 0.9))),
        !!other() := unname(quantile(tr[[other()]], c(0.5, 0.5, 0.9)))
      )
      plot_predictive_exceedance(model(), nd, cfg()$tail, lv)
    }
  })

  # ---- Trigger rain ----
  output$trigger_controls <- renderUI({
    req(!p1())
    ev <- fc(A()$event_probability)
    vals <- sort(unique(ev[[other()]]))
    sliderInput("other_val", paste("Condition on", var_label(other())),
      min = 0, max = signif(max(vals), 2), value = signif(vals[round(length(vals) / 4)], 2),
      step = signif(max(vals) / 50, 2), width = "100%")
  })
  output$trigger_left <- renderPlot({
    ev <- fc(A()$event_probability); tg <- fc(A()$triggers)
    req(ev, tg)
    prob <- as.numeric(input$prob)
    if (p1()) {
      plot_event_curve_1d(ev, tg, target(), prob) + labs(caption = cap())
    } else {
      req(input$other_val)
      vals <- sort(unique(ev[[other()]]))
      v <- vals[which.min(abs(vals - input$other_val))]
      plot_event_curve_1d(filter(ev, .data[[other()]] == v), filter(tg, .data[[other()]] == v), target(), prob) +
        labs(subtitle = sprintf("At %s = %.2f. Open circles: rain at which the chance reaches %s%%",
          var_label(other()), v, prob * 100), caption = cap())
    }
  })
  output$trigger_right <- renderPlot({
    if (p1()) {
      plot_predictive_exceedance(model(), tibble(!!target() := signif(seq(0.5, max(fc(A()$training)[[target()]]), length.out = 6), 2)),
        cfg()$tail, fc(A()$event_levels)) + labs(caption = cap())
    } else {
      plot_trigger_2d(fc(A()$triggers), target(), other(), rps = unique(c(2, rp(), 50))) +
        { if (!is.null(input$other_val)) geom_vline(xintercept = input$other_val, colour = ORANGE, linetype = "dashed") } +
        labs(caption = cap())
    }
  })
  output$trigger_table <- renderTable({
    tg <- fc(A()$triggers)
    req(tg)
    if (!p1()) {
      vals <- sort(unique(tg[[other()]]))
      pick <- vals[unique(round(seq(1, length(vals), length.out = 10)))]
      tg <- filter(tg, .data[[other()]] %in% pick, return_period == rp())
    }
    tg |>
      mutate(prob = paste0(prob * 100, "% chance"), trigger = ifelse(is.na(trigger), "not reached", sprintf("%.2f", trigger))) |>
      select(any_of(other()), `T (years)` = return_period, `Runoff level (mm/h)` = return_level, prob, trigger) |>
      pivot_wider(names_from = prob, values_from = trigger) |>
      rename_with(var_label, any_of(other()))
  }, digits = 2)

  # ---- Likeliest drivers ----
  output$likeliest_left <- renderPlot({
    lk <- fc(A()$likeliest)
    req(lk)
    if (p1()) plot_likeliest_1d(lk, target()) + labs(caption = cap())
    else plot_likeliest_2d_cond(lk, target(), other(), rp()) + labs(caption = cap())
  })
  output$likeliest_right <- renderPlot({
    lk <- fc(A()$likeliest)
    req(lk)
    if (p1()) {
      tr <- fc(A()$training)
      ggplot(tr, aes(.data[[target()]])) +
        geom_histogram(bins = 40, fill = "#86b6ef", colour = "white") +
        labs(x = var_label(target()), y = "Peaks", title = "Rainfall at observed runoff peaks", caption = cap()) +
        theme_ros()
    } else {
      d <- filter(lk, return_period == rp())
      ggplot(d, aes(.data[[target()]], .data[[other()]])) +
        geom_contour_filled(aes(z = f_given_event), bins = 8) +
        geom_point(data = fc(A()$training), colour = INK, alpha = 0.25, size = 0.6) +
        scale_fill_manual(values = grDevices::colorRampPalette(c("#f0efec", "#0d366b"))(8), guide = "none") +
        coord_cartesian(expand = FALSE) +
        labs(x = var_label(target()), y = var_label(other()),
          title = sprintf("Likeliest rain and snowmelt behind a %s-year peak", rp()),
          subtitle = "Darker = more likely, given the event occurred", caption = cap()) +
        theme_ros()
    }
  })

  # ---- Compare ----
  compare_data <- reactive({
    prob <- as.numeric(input$prob)
    purrr::map_dfr(list_analyses(root), function(nm) {
      p <- analysis_path(nm, "out", "triggers.rds", root = root)
      if (!file.exists(p)) return(NULL)
      c <- read_analysis(nm, root = root)
      tg <- readRDS(p) |> filter(cell_id == cell(), return_period == rp(), abs(.data$prob - .env$prob) < 1e-9)
      oth <- setdiff(names(tg), c("cell_id", "x", "y", "return_period", "return_level", "prob", "trigger"))
      tg |> mutate(analysis = nm, conditioning = if (length(oth)) .data[[oth[1]]] else NA_real_,
        conditioned_on = if (length(oth)) oth[1] else "(none)")
    })
  })
  output$compare_plot <- renderPlot({
    d <- compare_data()
    req(nrow(d) > 0)
    curves <- filter(d, !is.na(conditioning))
    flat <- filter(d, is.na(conditioning))
    ggplot() +
      geom_hline(data = flat, aes(yintercept = trigger, colour = analysis), linewidth = 1, linetype = "dashed") +
      geom_line(data = curves, aes(conditioning, trigger, colour = analysis), linewidth = 1) +
      scale_colour_manual(values = c("#2a78d6", ORANGE, "#1baf7a", "#4a3aa7"), name = NULL) +
      labs(x = "Snowmelt (mm/h), for analyses that condition on it", y = "Rain needed (mm/h)",
        title = sprintf("Rain needed for a %s%% chance of a %s-year runoff peak", as.numeric(input$prob) * 100, rp()),
        subtitle = "Dashed: analyses that ignore snowmelt. Solid: analyses conditioning on snowmelt.",
        caption = cap()) +
      theme_ros()
  })
  output$compare_table <- renderTable({
    compare_data() |>
      group_by(analysis, conditioned_on) |>
      summarise(`T-year level (mm/h)` = first(return_level),
        `Rain needed (no snowmelt / lowest)` = first(trigger[order(conditioning)]),
        `Rain needed (median snowmelt value)` = if (all(is.na(conditioning))) NA_real_ else
          trigger[which.min(abs(conditioning - median(unique(conditioning))))], .groups = "drop")
  }, digits = 2, na = "—")

  # ---- Diagnostics ----
  output$diag_rl <- renderPlot({ req(A()$return_levels); plot_return_levels(fc(A()$return_levels)) + facet_null() + labs(caption = cap()) })
  output$diag_pp <- renderPlot({ req(A()$diagnostics); plot_pp(A()$diagnostics) })
  output$diag_skill <- renderPlot({ req(A()$diagnostics); plot_skill(A()$diagnostics) })
}

shinyApp(ui, server)

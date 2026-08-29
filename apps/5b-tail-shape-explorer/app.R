# Tail shape explorer: is the shared tail index sensible, and how much does it
# matter?
#
# Run from the package root:
#   shiny::runApp("apps/5b-tail-shape-explorer")
#
# Requires:
#   - derived/era5_land_hourly_alps_dl_tail_summary.rds   (script 5)
#   - derived/era5_land_hourly_alps_dl_mixture_tails.rds  (script 5)
#
# Why this app exists. The cell marginal is an equal-weight mixture of per-hour
# predictive distributions, so its tail index is max_i(xi_i) unless the shape is
# shared across that cell's peak hours. Sharing it fixes that, but then the one
# shape per cell carries a lot of weight, and it is worth being able to see what
# it implies.
#
# EVERY PANEL IS WITHIN ONE CELL. The map is only there to pick a cell; nothing
# here compares a cell to its neighbours or borrows anything from them.
#
# Suggests: shiny, tidyverse, fs

library(shiny)
library(tidyverse)
library(fs)

repo_root <- here::here()
devtools::load_all(repo_root, quiet = TRUE)
source(file.path(repo_root, "apps", "ros_theme.R"))

summary_path <- path(repo_root, "derived", "era5_land_hourly_alps_dl_tail_summary.rds")
tails_path <- path(repo_root, "derived", "era5_land_hourly_alps_dl_mixture_tails.rds")

read_or_null <- function(p) if (file.exists(p)) readRDS(p) else NULL

tail_summary <- read_or_null(summary_path)
mixture_tails <- read_or_null(tails_path)

have_summary <- is.data.frame(tail_summary) && nrow(tail_summary) > 0
cells_ref <- if (have_summary) {
  distinct(tail_summary, cell_id, x, y)
} else {
  tibble(cell_id = integer(), x = double(), y = double())
}

tile_dims <- function(xy, x_col = "y", y_col = "x") {
  ux <- sort(unique(xy[[x_col]]))
  uy <- sort(unique(xy[[y_col]]))
  c(
    width = if (length(ux) > 1) stats::median(diff(ux)) else 0.25,
    height = if (length(uy) > 1) stats::median(diff(uy)) else 0.25
  )
}

nearest_cell <- function(lon, lat, tbl) {
  if (nrow(tbl) == 0) {
    return(NULL)
  }
  tbl$cell_id[which.min((tbl$y - lon)^2 + (tbl$x - lat)^2)]
}

REPORT_T <- rp_reporting()

ui <- fluidPage(
  theme = ros_bs_theme(),
  ros_header(
    "Tail shape",
    "The generalized-Pareto shape fitted to each cell's own peak hours, and what it costs to get wrong."
  ),
  uiOutput("banner"),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      selectInput(
        "cell_id",
        "Cell ID",
        choices = if (nrow(cells_ref)) sort(cells_ref$cell_id) else integer(),
        selected = if (nrow(cells_ref)) min(cells_ref$cell_id) else NULL
      ),
      radioButtons(
        "map_field",
        "Map shows",
        choices = c(
          "Tail shape (xi)" = "tail_shape",
          "Tail scale" = "tail_scale"
        ),
        selected = "tail_shape"
      ),
      selectInput(
        "sens_T",
        "Sensitivity: return period T (years)",
        choices = REPORT_T,
        selected = 200
      ),
      ros_note(
        "The sensitivity panel re-evaluates the selected cell's mixture over a ",
        "range of shapes, holding the scales fixed. The steeper that curve, the ",
        "more the design level rests on one estimated number."
      ),
      verbatimTextOutput("cell_meta")
    ),
    mainPanel(
      width = 9,
      fluidRow(
        column(7, plotOutput("map", height = "420px", click = "map_click")),
        column(5, plotOutput("return_curve", height = "420px"))
      ),
      fluidRow(
        column(6, plotOutput("sensitivity", height = "360px")),
        column(6, plotOutput("scale_spread", height = "360px"))
      )
    )
  )
)

server <- function(input, output, session) {
  output$banner <- renderUI({
    msgs <- character()
    if (!have_summary) {
      msgs <- c(msgs, paste(
        "Missing", basename(summary_path), "- run scripts/5-runoff_marginals.r."
      ))
    }
    if (!length(msgs)) {
      return(NULL)
    }
    div(class = "ros-note", lapply(msgs, tags$div))
  })

  observeEvent(input$map_click, {
    cid <- nearest_cell(input$map_click$x, input$map_click$y, cells_ref)
    if (!is.null(cid)) {
      updateSelectInput(session, "cell_id", selected = cid)
    }
  })

  selected_row <- reactive({
    req(have_summary)
    filter(tail_summary, cell_id == as.integer(input$cell_id))
  })

  selected_tail <- reactive({
    req(!is.null(mixture_tails))
    row <- filter(mixture_tails, cell_id == as.integer(input$cell_id))
    if (nrow(row) != 1L) {
      return(NULL)
    }
    list(mt = row$mixture_tail[[1]], nep = row$num_events_per_year[1])
  })

  output$map <- renderPlot({
    req(have_summary)
    td <- tile_dims(tail_summary)
    fld <- input$map_field
    d <- mutate(tail_summary, value = .data[[fld]], selected = cell_id == as.integer(input$cell_id))

    ggplot(d, aes(y, x, fill = value)) +
      geom_tile(width = td["width"], height = td["height"]) +
      geom_tile(
        data = filter(d, selected),
        fill = NA,
        colour = ros_palette$ink,
        linewidth = 0.7,
        width = td["width"],
        height = td["height"]
      ) +
      scale_fill_viridis_c(option = "F", direction = -1, na.value = "grey90") +
      coord_quickmap() +
      labs(
        title = if (fld == "tail_shape") {
          "Tail index xi by cell"
        } else {
          "Tail scale by cell"
        },
        subtitle = "Each cell fitted on its own data; click to select",
        fill = NULL
      ) +
      ros_ggtheme_map()
  })

  # The selected cell's own frequency-magnitude curve, straight off its mixture
  # tail. Nothing here involves any other cell.
  output$return_curve <- renderPlot({
    st <- selected_tail()
    if (is.null(st)) {
      return(
        ggplot() +
          annotate(
            "text",
            x = 0,
            y = 0,
            label = "No mixture tail for this cell",
            colour = ros_palette$muted
          ) +
          theme_void()
      )
    }
    rp <- exp(seq(log(2), log(500), length.out = 200))
    lvl <- mixture_tail_return_level(st$mt, rp * st$nep)
    d <- tibble(return_period = rp, level = lvl) |> filter(is.finite(level))
    validate(need(nrow(d) > 1, "This cell's tail region does not cover these return periods."))

    marks <- tibble(return_period = rp_reporting()) |>
      mutate(level = mixture_tail_return_level(st$mt, return_period * st$nep)) |>
      filter(is.finite(level))

    ggplot(d, aes(return_period, level)) +
      geom_line(colour = ros_palette$accent, linewidth = 0.9) +
      geom_point(data = marks, colour = ros_palette$warm, size = 2) +
      scale_x_log10() +
      labs(
        title = "This cell's runoff marginal",
        subtitle = "Closed-form mixture tail; dots are the reporting return periods",
        x = "Return period (years)",
        y = "Runoff (mm/h)"
      )
  })

  output$sensitivity <- renderPlot({
    st <- selected_tail()
    if (is.null(st)) {
      return(
        ggplot() +
          annotate("text", x = 0, y = 0, label = "No mixture tail for this cell",
                   colour = ros_palette$muted) +
          theme_void()
      )
    }
    mt <- st$mt
    tt <- as.numeric(input$sens_T)
    p <- 1 / (tt * st$nep)
    fitted_shape <- unique(round(mt$shape, 10))[1]

    grid <- seq(
      max(-0.3, fitted_shape - 0.35),
      min(0.9, fitted_shape + 0.35),
      length.out = 40
    )
    # Hold the scales fixed and sweep the shape, so the curve isolates the
    # leverage of the shape alone.
    lvl <- vapply(grid, function(xi) {
      alt <- mixture_tail(mt$threshold, mt$tail_prob, mt$scale, xi, mt$weights)
      mixture_tail_quantile(alt, p)
    }, numeric(1))

    at_fit <- mixture_tail_quantile(mt, p)
    ggplot(tibble(shape = grid, level = lvl), aes(shape, level)) +
      geom_line(colour = ros_palette$accent, linewidth = 0.9) +
      geom_vline(xintercept = fitted_shape, linetype = "dashed", colour = ros_palette$warm) +
      annotate(
        "point",
        x = fitted_shape,
        y = at_fit,
        colour = ros_palette$warm,
        size = 3
      ) +
      labs(
        title = sprintf("Leverage of the shape on the %s-year level", input$sens_T),
        subtitle = sprintf(
          "Scales held fixed; fitted shape %.3f gives %.1f mm/h",
          fitted_shape,
          at_fit
        ),
        x = "Tail shape xi",
        y = "Return level (mm/h)"
      )
  })

  # What the shared shape is riding on: how much this cell's own peak hours
  # differ in tail scale. A wide spread means the covariates are doing real work
  # within the cell; a narrow one means the hours look alike.
  output$scale_spread <- renderPlot({
    st <- selected_tail()
    row <- selected_row()
    if (is.null(st) || nrow(row) != 1L) {
      return(
        ggplot() +
          annotate(
            "text",
            x = 0,
            y = 0,
            label = "No mixture tail for this cell",
            colour = ros_palette$muted
          ) +
          theme_void()
      )
    }
    d <- tibble(
      scale = st$mt$scale,
      weight = st$mt$weights * st$mt$tail_prob
    ) |>
      filter(is.finite(scale), scale > 0, is.finite(weight), weight > 0)
    validate(need(nrow(d) > 1, "Only one tail component in this cell."))

    ggplot(d, aes(scale, weight = weight)) +
      geom_histogram(
        bins = 25,
        fill = ros_palette$accent_soft,
        colour = ros_palette$accent,
        linewidth = 0.25
      ) +
      geom_vline(
        xintercept = row$tail_scale,
        colour = ros_palette$warm,
        linewidth = 1
      ) +
      labs(
        title = "Spread of tail scale within this cell",
        subtitle = sprintf(
          "%d components; the line is the single GP scale that matches the cell in the limit (%.3g)",
          nrow(d),
          row$tail_scale
        ),
        x = "GP scale of a component",
        y = "Weight"
      )
  })

  output$cell_meta <- renderText({
    row <- selected_row()
    if (nrow(row) != 1L) {
      return("No cell selected.")
    }
    st <- selected_tail()
    lines <- c(
      sprintf("cell_id       %d", row$cell_id),
      sprintf("lon, lat      %.3f, %.3f", row$y, row$x),
      sprintf("events / year %.2f", row$num_events_per_year),
      sprintf("tail shape    %.4f", row$tail_shape),
      sprintf("tail scale    %.3f", row$tail_scale),
      sprintf("threshold     %.3f", row$tail_threshold)
    )
    if (!is.null(st)) {
      lines <- c(
        lines,
        sprintf("components    %d", length(st$mt$threshold)),
        sprintf(
          "valid above   %.2f mm/h (p <= %.3f)",
          mixture_tail_lower(st$mt),
          mixture_tail_max_prob(st$mt)
        )
      )
    }
    paste(lines, collapse = "\n")
  })
}

shinyApp(ui, server)

# Tail shape explorer: is the shared tail index sensible, and how much does it
# matter?
#
# Run from the package root:
#   shiny::runApp("apps/5b-tail-shape-explorer")
#
# Requires (all from scripts/4 and 5):
#   - derived/era5_land_hourly_alps_dl_tail_summary.rds   (script 5)
#   - derived/era5_land_hourly_alps_dl_mixture_tails.rds  (script 5)
#   - derived/era5_land_hourly_alps_dl_tail_shapes.rds    (script 4, when
#     spatial smoothing is on -- optional)
#
# Why this app exists. The cell marginal is an equal-weight mixture of per-hour
# predictive distributions, so its tail index is max_i(xi_i) unless the shape is
# shared. Sharing it fixes that, but then the single shape per cell carries a
# lot of weight, and it is worth being able to see it: whether it varies
# smoothly over the map, how far smoothing moved it, and how much the design
# return level actually depends on it.
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
shapes_path <- path(repo_root, "derived", "era5_land_hourly_alps_dl_tail_shapes.rds")

read_or_null <- function(p) if (file.exists(p)) readRDS(p) else NULL

tail_summary <- read_or_null(summary_path)
mixture_tails <- read_or_null(tails_path)
shape_pairs <- read_or_null(shapes_path)

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
    "The shared generalized-Pareto shape per cell: where it varies, how much smoothing moved it, and what it costs to get wrong."
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
        column(5, plotOutput("shrinkage", height = "420px"))
      ),
      fluidRow(
        column(6, plotOutput("sensitivity", height = "360px")),
        column(6, plotOutput("neighbours", height = "360px"))
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
    if (is.null(shape_pairs)) {
      msgs <- c(msgs, paste(
        "No", basename(shapes_path),
        "- spatial smoothing is off in inputs/distributional_learning.yaml,",
        "so the raw-vs-smoothed panel is empty."
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
        subtitle = "Click a cell to select it",
        fill = NULL
      ) +
      ros_ggtheme_map()
  })

  output$shrinkage <- renderPlot({
    if (is.null(shape_pairs)) {
      return(
        ggplot() +
          annotate(
            "text",
            x = 0,
            y = 0,
            label = "Spatial smoothing not enabled",
            colour = ros_palette$muted
          ) +
          theme_void()
      )
    }
    sel <- as.integer(input$cell_id)
    d <- mutate(shape_pairs, selected = cell_id == sel)
    ggplot(d, aes(shape_raw, shape_smoothed)) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = ros_palette$muted) +
      geom_point(colour = ros_palette$accent, alpha = 0.6, size = 1.8) +
      geom_point(
        data = filter(d, selected),
        colour = ros_palette$warm,
        size = 3.4
      ) +
      labs(
        title = "How far smoothing moved each cell",
        subtitle = "Points on the dashed line kept their own estimate",
        x = "Shape before smoothing",
        y = "Shape after smoothing"
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

  output$neighbours <- renderPlot({
    req(have_summary)
    sel <- as.integer(input$cell_id)
    row <- selected_row()
    req(nrow(row) == 1L)

    nb <- grid_neighbours(tail_summary$x, tail_summary$y, radius = 1)
    i <- which(tail_summary$cell_id == sel)
    idx <- nb[[i]]

    d <- tibble(
      what = c("this cell", rep("neighbours", length(idx))),
      shape = c(row$tail_shape, tail_summary$tail_shape[idx])
    )
    ggplot(d, aes(what, shape)) +
      geom_hline(
        yintercept = stats::median(tail_summary$tail_shape, na.rm = TRUE),
        linetype = "dotted",
        colour = ros_palette$muted
      ) +
      geom_point(
        aes(colour = what),
        size = 3,
        alpha = 0.8,
        position = position_jitter(width = 0.08, height = 0)
      ) +
      scale_colour_manual(
        values = c("this cell" = ros_palette$warm, "neighbours" = ros_palette$accent)
      ) +
      labs(
        title = "This cell against its neighbours",
        subtitle = "Dotted line is the grid-wide median shape",
        x = NULL,
        y = "Tail shape xi"
      ) +
      theme(legend.position = "none")
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
